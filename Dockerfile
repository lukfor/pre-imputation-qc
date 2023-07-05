FROM nfcore/base:1.14
MAINTAINER Lukas Forer <lukas.forer@i-med.ac.at>

COPY environment.yml .
RUN conda env create --quiet -f environment.yml && conda clean -a

ENV PATH /opt/conda/envs/pre-imputation/bin:$PATH


# Install jbang (not as conda package available)
WORKDIR "/opt"
RUN wget https://github.com/jbangdev/jbang/releases/download/v0.59.0/jbang.zip && \
    unzip -q jbang.zip && \
    rm jbang.zip
ENV PATH="/opt/jbang/bin:${PATH}"

# Install genomic-utils
WORKDIR "/opt"
ENV GENOMIC_UTILS_VERSION="v0.2.0"
RUN wget https://github.com/genepi/genomic-utils/releases/download/${GENOMIC_UTILS_VERSION}/genomic-utils.jar

# Install imputation bot
ENV IMPUTATIONBOT_VERSION="0.9.4"
RUN mkdir /opt/imputationbot
WORKDIR "/opt/imputationbot"
RUN wget https://github.com/lukfor/imputationbot/releases/download/v${IMPUTATIONBOT_VERSION}/imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    unzip -q imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    rm imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    ./imputationbot version
ENV PATH="/opt/imputationbot:${PATH}"
