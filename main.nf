#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_sample_list } from './modules/irida-uploader.nf'
include { irida_uploader } from './modules/irida-uploader.nf'

workflow {

  ch_config = Channel.fromPath("${params.config}")
  ch_run_dir = Channel.fromPath("${params.run_dir}", type: 'dir')

  if(params.instrument_type == "nextseq") {
    fastq_subdir = "Analysis/1/Data/fastq"
  } else {
    fastq_subdir = "Data/Intensities/BaseCalls"
  }
    
  ch_fastq = Channel.fromFilePairs("${params.run_dir}/${fastq_subdir}/*_{R1,R2}_*.fastq.gz")

  main:
    create_sample_list(ch_run_dir)
    ch_fastq_flattened = ch_fastq.map{ it -> [ it[1][0], it[1][1] ]}.flatten()
    // irida_uploader(create_sample_list.out.combine(ch_config), ch_fastq_flattened.collect())
}
