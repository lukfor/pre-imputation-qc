# Pre-Imputation Pipeline

> Nextflow pipeline to prepare your data for Genotype Imputation

## Requirements

Before running this pipeline, make sure you have the following dependencies installed:

- [Nextflow](https://www.nextflow.io/)

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

Create a config file (e.g. `my-project.config` in your project folder) and set the paths to your data:

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
nextflow run main.nf -c my-project.config
```

A html report is created in `output` and the `vcf.gz` are ready to submit them to the Michigan Imputation Server.

## Parameters

Here's a table of the parameters along with their descriptions and default values:

| Parameter         | Description                                  | Default Value      | Required |
| ----------------- | -------------------------------------------- | ------------------ | -------- |
| project           | Name of the project                          | null               | Yes      |
| input             | Directory containing input PLINK files       | null               | Yes      |
| output            | Output directory for pipeline results        | "output/genotyped" | Yes      |
| chip              | Type of genotyping chip used                 | null               | Yes      |
| build             | Reference genome build (e.g., hg19, hg38)    | null               | Yes      |
| strand_file       | Path to the strand file                      | null               | Yes      |
| refalt_file       | Path to the reference/alternate allele file  | null               | Yes      |
| chunkSize         | Chunk size for processing data               | 20000000           | No       |
| minSampleCallRate | Minimum sample call rate threshold           | 0.5                | No       |
| minSnpCallRate    | Minimum SNP call rate threshold              | 0.9                | No       |
| maf               | Minor allele frequency threshold             | 0                  | No       |
| hwe               | Hardy-Weinberg equilibrium p-value threshold | 1E-6               | No       |
| cleanSampleIds    | Clean sample IDs (true/false)                | false              | No       |
| excludeSamples    | File containing samples to exclude           | null               | No       |
| useDoubleId       | Use double IDs for samples (true/false)      | true               | No       |

These parameters allow customization of the pipeline's behavior and thresholds for various data processing steps. You can adjust these values as needed to fit your specific project requirements.

## Workflow

1. **Clean Sample IDs**: (Optional) Clean sample IDs in input files if specified.

2. **Exclude Samples**: (Optional) Exclude specified samples from the analysis.

3. **Filter and Fix Strand Flips**: Filter and fix strand flips using provided strand and ref/alt allele files.

4. **Merge VCF Files**: Merge filtered VCF files.

5. **Filter Merged VCF**: Further filter the merged VCF file.

6. **Create Final PLINK**: Generate the final PLINK files.

7. **Split Into Chromosomes**: Split the final PLINK files into separate chromosomes.

8. **Create Report**: Generate a report summarizing the pipeline's execution.

## Outputs

- Processed PLINK files for each chromosome.
- Various reports and statistics files.

## Contact

Lukas Forer (@lukfor), Institute of Genetic Epidemiology, Medical University of Innsbruck

## License

This pipeline is distributed under the [MIT License](LICENSE).
