FROM continuumio/miniconda
MAINTAINER Lukas Forer <lukas.forer@i-med.ac.at>

COPY environment.yml .
RUN \
   conda env update -n root -f environment.yml \
&& conda clean -a

# Install jbang (not as conda package available)
WORKDIR "/opt"
RUN wget https://github.com/jbangdev/jbang/releases/download/v0.59.0/jbang.zip && \
    unzip -q jbang.zip && \
    rm jbang.zip
ENV PATH="/opt/jbang/bin:${PATH}"

# Install imputation bot
ENV IMPUTATIONBOT_VERSION="0.9.4"
RUN mkdir /opt/imputationbot
WORKDIR "/opt/imputationbot"
RUN wget https://github.com/lukfor/imputationbot/releases/download/v${IMPUTATIONBOT_VERSION}/imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    unzip -q imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    rm imputationbot-${IMPUTATIONBOT_VERSION}-linux.zip && \
    ./imputationbot version
ENV PATH="/opt/imputationbot:${PATH}"
