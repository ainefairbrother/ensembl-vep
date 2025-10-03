#!/usr/bin/env nextflow

/* 
 * Script to check if the files are bgzipped and bgzip if not
 */

import java.io.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

def checkVCFheader (f) {
  // Check file extension
  if (!(f ==~ /.*\.vcf$/) && !(f ==~ /.*\.vcf\.b?gz$/)) {
    return false
  }

  // Read lines (compressed or not)
  List<String> lines
  if (f ==~ /.*\.b?gz$/) {
    def br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))))
    lines = []
    for (def L = br.readLine(); L != null; L = br.readLine()) { lines << L }
    br.close()
  } else {
    lines = new File(f).readLines()
  }

  // Check file header
  def is_vcf_format = false
  def has_header = false
  for (String line : lines) {
    if (!line.startsWith('#')) {
      // stop inspecting file when reaching a line not starting with hash
      break
    } else if (line.startsWith('##fileformat=')) {
      is_vcf_format = true
    } else if (line.startsWith('#CHROM')) {
      has_header = true
    }
  }
  return is_vcf_format && has_header
}

process checkVCF {
  /*
  Function to check input VCF files

  Returns
  -------
  Tuple of VCF, VCF index, vep config file, a output dir, and the index type of VCF file
  */

  cpus params.cpus
  label 'vep'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(vcf), path(vcf_index), path(vep_config)
  
  output:
  tuple val(meta), path("*.{gz,bgz}", includeInputs: true), path ("*.{gz,bgz}.{tbi,csi}", includeInputs: true), path(vep_config)

  afterScript "rm -f *.vcf *.vcf.tbi *.vcf.csi tmp.vcf"

  script:
  def index_type = meta.index_type
  def tabix_arg  = index_type == 'tbi' ? '' : '-C'
  def isGzipped  = (vcf.extension in ['gz','bgz'])
  def out_vcf    = isGzipped ? "${vcf}" : "${vcf}.gz"

  def sort_cmd = ""
  if (params.sort) {
    def cat_cmd = isGzipped ? "zcat ${vcf}" : "cat ${vcf}"
    sort_cmd  = "(${cat_cmd} | head -1000 | grep '^#'; ${cat_cmd} | grep -v '^#' | sort -k1,1d -k2,2n) > tmp.vcf; "
    sort_cmd += isGzipped ? "bgzip -c tmp.vcf > ${vcf}; rm -f tmp.vcf;" : "mv -f tmp.vcf ${vcf};"
  }

  """
  ${sort_cmd}

  # Compress only if not already gz/bgz
  if [ "${isGzipped}" != "true" ]; then
    bgzip -c ${vcf} > ${out_vcf}
  fi

  # fix accidental .gz.gz -> .gz (harmless if none)
  for f in ./*.gz.gz; do
    [ -e "\$f" ] || break
    mv -vn -- "\$f" "\${f%.gz}"
  done

  # Index the specific output file (not a glob)
  [ -f "${out_vcf}.${index_type}" ] || tabix ${tabix_arg} -p vcf -f "${out_vcf}"
  """
}