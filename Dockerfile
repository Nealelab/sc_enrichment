FROM debian:jessie
ENV CLOUD_SDK_VERSION 193.0.0

MAINTAINER Andrea Ganna

ARG INSTALL_COMPONENTS
RUN apt-get update -qqy && apt-get install -qqy \
        curl \
        gcc \
        python-dev \
        python-setuptools \
        apt-transport-https \
        lsb-release \
        openssh-client \
        git \
    && easy_install -U pip && \
    pip install -U crcmod && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version

RUN apt-get install -y unzip wget bedtools zlib1g-dev g++

RUN pip install -U joblib==0.11 pandas==0.19.2 numpy==1.11.3 scipy==0.18.1 bitarray==0.8.1 pybedtools==0.7.10 

RUN	mkdir -p /home/mtag/ && \
    wget --quiet -P /home/mtag/ https://github.com/omeed-maghzian/mtag/archive/master.zip && \
    unzip -q /home/mtag/master.zip -d /home/mtag

RUN	mkdir -p /home/ldscore/ && \
	wget --quiet -P /home/ldscore/ https://github.com/bulik/ldsc/archive/master.zip && \
    unzip -q /home/ldscore/master.zip -d /home/ldscore

RUN	mkdir -p /home/sc_enrichement/ && \
	wget --quiet -P /home/sc_enrichement/ https://github.com/Nealelab/sc_enrichement/archive/master.zip && \
    unzip -q /home/sc_enrichement/master.zip -d /home/sc_enrichement

VOLUME ["/root/.config"]

CMD [ "/bin/bash" ]