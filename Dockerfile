FROM python:2.7-slim

MAINTAINER Andrea Ganna


ENV GOOGLE_PROJECT ukbb-gay-was
ENV PATH $PATH:/root/google-cloud-sdk/bin

RUN apt-get update --fix-missing && \
	apt-get install -y unzip python-pip wget curl bedtools zlib1g-dev git

RUN pip install google-compute-engine joblib==0.11 pandas==0.19.2 numpy==1.11.3 scipy==0.18.1 bitarray==0.8.1 pybedtools==0.7.10


RUN	mkdir -p /home/mtag/ && \
    wget --quiet -P /home/mtag/ https://github.com/omeed-maghzian/mtag/archive/master.zip && \
    unzip -q /home/mtag/master.zip -d /home/mtag

RUN	mkdir -p /home/ldscore/ && \
	wget --quiet -P /home/ldscore/ https://github.com/bulik/ldsc/archive/master.zip && \
    unzip -q /home/ldscore/master.zip -d /home/ldscore


RUN	mkdir -p /home/sc_enrichement/ && \
	wget --quiet -P /home/sc_enrichement/ https://github.com/Nealelab/sc_enrichement/archive/master.zip && \
    unzip -q /home/sc_enrichement/master.zip -d /home/sc_enrichement


RUN curl -sSL https://sdk.cloud.google.com | bash

ADD UKBB_GAY_WAS-d00c97def1c8.json /root/service_key.json

RUN gcloud auth activate-service-account --key-file /root/service_key.json

RUN gcloud config set compute/zone us-central1-f

RUN gcloud config set container/cluster ukbb-gay-was

CMD [ "/bin/bash" ]
