#!/bin/bash

VERSION="v1.2.0"
IMAGE="genepi/pre-imputation-qc"
IMAGE_FILENAME="genepi-pre-imputation-qc"

mkdir -p containers

docker build -t ${IMAGE}:${VERSION} . --platform linux/amd64
#docker save -o containers/${IMAGE_FILENAME}-${VERSION}.docker ${IMAGE}:${VERSION} 
singularity build containers/${IMAGE_FILENAME}-${VERSION}.sif docker-daemon://${IMAGE}:${VERSION}
