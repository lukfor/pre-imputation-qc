# Pre-Imputation Pipeline

> Nextflow pipeline to prepare your data for Genotype Imputation

## Requirements

- Nextflow:

```
curl -s https://get.nextflow.io | bash
```

- Docker

## Installation

Build docker image before run the pipeline:

```
docker build -t lukfor/pre-imputation . # don't ingore the dot here
```


Test the pipeline and the created docker image with test-data:

```
nextflow run main.nf
```

## Usage

Create a config file (e.g. `pre-imputation-qc.config` in your project folder) and set the paths to your data:

```
params.project = "my-project"
params.input = "data/*/*.{map,ped}"
params.output = "output"
params.chip = "GSAMD-24v3-0-EA_20034606_A1.b37"
params.strand_file = "data/${params.chip}.strand"
params.refalt_file = "data/${params.chip}.RefAlt"

params.chunkSize= 20000000
params.minSampleCallRate = 0.5
params.minSnpCallRate = 0.9
```

Execute it the pipeline with you project specific configuration from your project folder:

```
nextflow run /home/lukas/git/pre-imputation-qc/main.nf -c pre-imputation-qc.config
```

A html report is created in `output` and the `vcf.gz` are ready to submit them to the Michigan Imputation Server.

## Contact

Lukas Forer (@lukfor), Institute of Genetic Epidemiology, Medical University of Innsbruck

## License
gwas-scripts is MIT Licensed.
