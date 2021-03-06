# Routine IRIDA Upload

Pipeline for performing routine uploads of illumina sequence data to the [IRIDA](http://www.irida.ca) analysis platform.

## Usage
```
nextflow run BCCDC-PHL/routine-irida-upload --run_dir <illumina run directory> --outdir <output directory> [--config <irida uploader config>]
```

## Outputs
```
<outdir>
├── <run_id>_irida_uploader.log
└── <run_id>_irida_uploader_status.info
```