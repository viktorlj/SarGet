#!/bin/sh
set -eu pipefail

# Local ENV variables
BCFTOOLS_VERSION=1.5
BWA_VERSION=0.7.16
HTSLIB_VERSION=1.5
SAMTOOLS_VERSION=1.5

# Install libraries
apk update && apk add \
  alpine-sdk \
  bzip2 \
  ca-certificates \
  curl-dev \
  g++ \
  gcc \
  git \
  gnupg \
  bzip2-dev \
  xz-dev \
  ncurses-dev \
  openssl-dev \
  make \
  python \
  python3 \
  unzip \
  wget \
  zlib-dev

# Install tools
mkdir /build

# Install BBduk
cd /build
wget --quiet https://kent.dl.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz
tar xvfz BBMap_38.44.tar.gz
mkdir /bbmap
mv bbmap/* /bbmap

# Install fgbio
cd /build
wget --quiet https://github.com/fulcrumgenomics/fgbio/releases/download/0.8.1/fgbio-0.8.1.jar
mkdir /fgbio
mv fgbio-0.8.1.jar /fgbio/

# Install BWA
cd /build
git clone http://github.com/lh3/bwa.git bwa \
  --branch v${BWA_VERSION}
cd bwa
make
cp bwa /usr/local/bin/bwa

# Install HTSlib
cd /build
wget --quiet -O htslib-${HTSLIB_VERSION}.tar.bz2 \
  https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
tar xfj htslib-${HTSLIB_VERSION}.tar.bz2
rm htslib-${HTSLIB_VERSION}.tar.bz2
cd htslib-${HTSLIB_VERSION}
./configure --prefix=/opt/htslib
make && make install

# Install pindel
cd /build
git clone https://github.com/genome/pindel.git
cd pindel
./INSTALL /opt/htslib/
mkdir /pindel
mv /build/pindel/* /pindel/

# Install BCFtools
cd /build
wget --quiet -O bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
  https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar xfj bcftools-${BCFTOOLS_VERSION}.tar.bz2
rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd bcftools-${BCFTOOLS_VERSION}
./configure --prefix=/opt/bcftools
make && make install

# Install Samtools
cd /build
wget --quiet -O samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xfj samtools-${SAMTOOLS_VERSION}.tar.bz2
rm samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}
./configure --prefix=/opt/samtools
make && make install

#Install primerclip (specific for Swift data)
cd /build
git clone https://github.com/swiftbiosciences/primerclip primerclip
cd primerclip
mkdir /primerclip
cp .stack-work/install/x86_64-linux/lts-11.0/8.2.2/bin/primerclip /primerclip/

#Install Agent toolkit
cd /build
wget --quiet -O AGeNT.zip https://dt4ei3l3hxs7z.cloudfront.net/
unzip AGeNT.zip
rm AGeNT.zip
mkdir /AGeNT
mv LocatIt_v4.0.1.jar /AGeNT/LocatIt.jar
mv SurecallTrimmer_v4.0.1.jar /AGeNT/SurecallTrimmer.jar

# Install lofreq
cd /build
wget --quiet https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.3.1_linux-x86-64.tgz
tar xvfz lofreq_star-2.1.3.1_linux-x86-64.tgz
mkdir /lofreq
mv lofreq_star-2.1.3.1/bin/* /lofreq

# Clean up install
cd /

rm -rf /build /var/lib/apt/lists/* /opt/get-pip.py