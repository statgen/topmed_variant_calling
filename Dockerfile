FROM ubuntu:16.04

COPY . /topmed_variant_calling

RUN apt-get update && apt-get install -y \
  apt-utils \
  automake \
  autoconf \
  build-essential \
  git \
  ghostscript \
  gnuplot \
  groff \
  libcurl4-openssl-dev \
  liblzma-dev \
  libncurses5-dev \
  libssl-dev \
  libzstd-dev \
  r-base \
  unzip \
  wget \
  zlib1g-dev

RUN mkdir /tmp/plink && cd /tmp/plink && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && unzip plink_linux_x86_64_20190617.zip && cp plink /usr/local/bin/plink-1.9

WORKDIR /topmed_variant_calling
RUN rm -r /tmp/plink

RUN git submodule init && git submodule update 

RUN cd libsvm/ && make clean && make && cd ..
RUN cd apigenome && autoreconf -vfi && ./configure --prefix $PWD && make clean && make && make install && cd ..
RUN cd libStatGen && make clean && make && cd ..
RUN cd bamUtil && make clean && make && cd ..
RUN cd invNorm && make clean && make && cd ..
RUN cd htslib && autoheader && autoconf && ./configure && make clean && make && cd ..
RUN cd vt-topmed && make clean && make && cd ..
RUN cd cramore && autoreconf -vfi && ./configure && make clean && make && cd ..
RUN cd samtools && autoheader && autoconf -Wno-syntax && ./configure && make clean && make && cd ..
RUN cd bcftools && make clean && make && cd ..
RUN cd king && rm -f *.o && g++ -O3 -c *.cpp && g++ -O3 -o king *.o -lz && cd ..

