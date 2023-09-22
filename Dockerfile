FROM nfcore/base:1.14
MAINTAINER Lukas Forer <lukas.forer@i-med.ac.at>

COPY environment.yml .
RUN conda env create --quiet -f environment.yml && conda clean -a

ENV PATH /opt/conda/envs/pre-imputation/bin:$PATH

# Install genomic-utils
WORKDIR "/opt"
ENV GENOMIC_UTILS_VERSION="v0.2.0"
RUN wget https://github.com/genepi/genomic-utils/releases/download/${GENOMIC_UTILS_VERSION}/genomic-utils.jar
