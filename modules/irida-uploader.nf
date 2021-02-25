process create_sample_list {

  tag { run_id }

  executor 'local'

  input:
  path(run_dir)

  output:
  tuple val(run_id), path("SampleList.csv")

  script:
  run_id = run_dir.baseName
  """
  cp ${run_dir}/SampleSheet.csv .
  echo "[Data]" > SampleList_Header.csv
  echo "Sample_Name,Project_ID,File_Forward,File_Reverse" >> SampleList_Header.csv
  sed -e '1,/\\[Data\\]/d' SampleSheet.csv > SampleSheet_Data.csv
  tail -n+2 SampleSheet_Data.csv | awk -F "," 'BEGIN { OFS=FS }; \$10 { print \$1, \$10 }' > Sample_Name_Project_ID.csv
  while IFS="," read -r sample_name project_id; do \
    ls -1 ${run_dir}/Data/Intensities/BaseCalls/\${sample_name}*R1*.fastq.gz | xargs -n 1 basename >> Reads_R1.csv; \
  done < Sample_Name_Project_ID.csv
  while IFS="," read -r sample_name project_id; do \
    ls -1 ${run_dir}/Data/Intensities/BaseCalls/\${sample_name}*R2*.fastq.gz | xargs -n 1 basename >> Reads_R2.csv; \
  done < Sample_Name_Project_ID.csv
  paste -d "," Sample_Name_Project_ID.csv Reads_R1.csv Reads_R2.csv > SampleList_Body.csv
  cat SampleList_Header.csv SampleList_Body.csv > SampleList.csv
  """
}

process irida_uploader {

  tag { run_id }

  executor 'local'

  input:
  tuple val(run_id), path(sample_list), path(config)
  path(reads)

  output:
  tuple val(run_id), path("${run_id}_irida_uploader.log"), path("${run_id}_irida_uploader_status.info")

  script:
  def config = config.name != 'NO_FILE' ? "-c ${config}" : ""
  """
  irida-uploader ${config} --config_parser directory -d .
  mv irida-uploader.log ${run_id}_irida_uploader.log
  mv irida_uploader_status.info ${run_id}_irida_uploader_status.info
  """
}
