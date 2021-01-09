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


## Usage

Test pipeline with test-data:

```
nextflow run main.nf
```

## Contact

Lukas Forer (@lukfor), Institute of Genetic Epidemiology, Medical University of Innsbruck

## License
gwas-scripts is MIT Licensed.
