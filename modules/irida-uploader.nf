process create_sample_list {

  tag { run_id }

  executor 'local'

  input:
  path(run_dir)

  output:
  tuple val(run_id), path("SampleList.csv")

  script:
  run_id = run_dir.baseName
  if (params.instrument_type == "nextseq") {
    data_prefix = "Cloud_"
    awk_string = "BEGIN { OFS=FS }; \$2 ~ /^[[:digit:]]+\$/ { print \$1, \$2 }"
  } else {
    data_prefix = ""
    awk_string = "BEGIN { OFS=FS }; \$10 ~ /^[[:digit:]]+\$/ { print \$2, \$10 }"
  }
  """
  cp ${run_dir}/SampleSheet*.csv SampleSheet_Copy.csv
  echo "[Data]" > SampleList_Header.csv
  echo "Sample_Name,Project_ID,File_Forward,File_Reverse" >> SampleList_Header.csv
  sed -e '1,/\\[${data_prefix}Data\\]/d' SampleSheet_Copy.csv > SampleSheet_Data.csv
  tail -n+2 SampleSheet_Data.csv | awk -F "," '${awk_string}' > Sample_Name_Project_ID.csv
  touch Reads_R1.csv
  while IFS="," read -r sample_name project_id; do \
    ls -1 ${run_dir}/${params.fastq_subdir}/\${sample_name}*R1*.fastq.gz | xargs -n 1 basename >> Reads_R1.csv; \
  done < Sample_Name_Project_ID.csv
  touch Reads_R2.csv
  while IFS="," read -r sample_name project_id; do \
    ls -1 ${run_dir}/${params.fastq_subdir}/\${sample_name}*R2*.fastq.gz | xargs -n 1 basename >> Reads_R2.csv; \
  done < Sample_Name_Project_ID.csv
  paste -d "," Sample_Name_Project_ID.csv Reads_R1.csv Reads_R2.csv > SampleList_Body.csv
  cat SampleList_Header.csv SampleList_Body.csv > SampleList.csv
  """
}

process irida_uploader {

  tag { run_id }

  executor 'local'

  publishDir "${params.outdir}", pattern: "${run_id}_irida_uploader*", mode: 'copy'
  publishDir "${params.outdir}", pattern: "pipeline_complete.json", mode: 'copy'

  input:
    tuple val(run_id), path(sample_list), path(config)
    path(reads)

  output:
    tuple val(run_id), path("${run_id}_irida_uploader.log"), path("${run_id}_irida_uploader_status.info"), path("pipeline_complete.json")

  script:
    def config = config.name != 'NO_FILE' ? "-c ${config}" : ""
    """
    irida-uploader ${config} --config_parser directory -d . || true
    mv irida-uploader.log ${run_id}_irida_uploader.log
    mv irida_uploader_status.info ${run_id}_irida_uploader_status.info
    echo "{\\"pipeline_name\\": \\"BCCDC-PHL/routine-irida-upload\\", \\"timestamp_completed\\": \\"\$(date --iso=seconds)\\"}" > pipeline_complete.json
    """
}
