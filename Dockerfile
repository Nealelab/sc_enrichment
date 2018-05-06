FROM google/cloud-sdk:slim

MAINTAINER Andrea Ganna

RUN apt-get install -y unzip wget bedtools zlib1g-dev g++

RUN pip install -U joblib==0.11 pandas==0.19.2 numpy==1.11.3 scipy==0.18.1 bitarray==0.8.1 pybedtools==0.7.10 h5py==2.7.1

RUN	mkdir -p /home/mtag/ && \
    wget --quiet -P /home/mtag/ https://github.com/omeed-maghzian/mtag/archive/master.zip && \
    unzip -q /home/mtag/master.zip -d /home/mtag

RUN	mkdir -p /home/ldscore/ && \
	wget --quiet -P /home/ldscore/ https://github.com/Nealelab/ldsc/archive/kt_exclude_files.zip && \
    unzip -q /home/ldscore/kt_exclude_files.zip -d /home/ldscore

RUN	mkdir -p /home/sc_enrichement/ && \
	wget --quiet -P /home/sc_enrichement/ https://github.com/Nealelab/sc_enrichement/archive/master.zip && \
    unzip -q /home/sc_enrichement/master.zip -d /home/sc_enrichement

VOLUME ["/root/.config"]

CMD [ "/bin/bash" ]