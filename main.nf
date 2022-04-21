#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_sample_list } from './modules/irida-uploader.nf'
include { irida_uploader } from './modules/irida-uploader.nf'

workflow {

  ch_config = Channel.fromPath("${params.config}")
  ch_run_dir = Channel.fromPath("${params.run_dir}", type: 'dir')

  ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }


  main:
    create_sample_list(ch_run_dir)
    ch_fastq_flattened = ch_fastq.map{ it -> [ it[1], it[2] ]}.flatten()
    irida_uploader(create_sample_list.out.combine(ch_config), ch_fastq_flattened.collect())
}
