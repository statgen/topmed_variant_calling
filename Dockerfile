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
  python3 \
  r-base \
  unzip \
  wget \
  zlib1g-dev

RUN mkdir /tmp/plink && cd /tmp/plink && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && unzip plink_linux_x86_64_20190617.zip && cp plink /usr/local/bin/plink-1.9

WORKDIR /topmed_variant_calling
RUN rm -r /tmp/plink

RUN git submodule init && git submodule update 

RUN cd libsvm/ && git clean -fdx && make && cd ..
RUN cd apigenome && git clean -fdx && autoreconf -vfi && ./configure --prefix $PWD && make && make install && cd ..
RUN cd libStatGen && git clean -fdx && make && cd ..
RUN cd bamUtil && git clean -fdx && make && cd ..
RUN cd invNorm && git clean -fdx && make && cd ..
RUN cd htslib && git clean -fdx && autoheader && autoconf && ./configure && make && cd ..
RUN cd vt-topmed && git clean -fdx && make && cd ..
RUN cd cramore && git clean -fdx && autoreconf -vfi && ./configure && make && cd ..
RUN cd samtools && git clean -fdx && autoheader && autoconf -Wno-syntax && ./configure && make && cd ..
RUN cd bcftools && git clean -fdx && make && cd ..
RUN cd king && rm -f king *.o && g++ -O3 -c *.cpp && g++ -O3 -o king *.o -lz && cd ..

