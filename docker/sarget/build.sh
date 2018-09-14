#!/bin/sh
set -eu pipefail

# Local ENV variables
BCFTOOLS_VERSION=1.5
BWA_VERSION=0.7.16
HTSLIB_VERSION=1.5
SAMTOOLS_VERSION=1.5
SAMBAMBA_VERSION=0.6.7

# Install libraries
apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  bzip2 \
  ca-certificates \
  g++ \
  gcc \
  git \
  libbz2-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libssl-dev \
  make \
  python \
  python3 \
  unzip \
  wget \
  zlib1g-dev

# Install dotnet
wget -qO- https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor > microsoft.asc.gpg
mv microsoft.asc.gpg /etc/apt/trusted.gpg.d/
wget -q https://packages.microsoft.com/config/debian/8/prod.list
mv prod.list /etc/apt/sources.list.d/microsoft-prod.list
apt-get install -y --no-install-recommends apt-transport-https
apt-get update
apt-get install -y --no-install-recommends dotnet-sdk-2.1

# Install tools
mkdir /build

# Install sambamba
cd /build
wget --quiet -O sambamba-${SAMBAMBA_VERSION}.tar.bz2 \
  https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2
tar xfj sambamba-${SAMBAMBA_VERSION}.tar.bz2
cp sambamba /usr/local/bin/sambamba

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

#Install libStatGen
cd /build
git clone https://github.com/statgen/libStatGen.git libStatGen
cd libStatGen
make

#Install bam-utils
cd /build
git clone https://github.com/statgen/bamUtil.git bamUtil
cd bamUtil
make all
cp bin/bam /usr/local/bin/bam

#Install Agent toolkit
cd /build
wget --quiet -O AGeNT.zip https://dt4ei3l3hxs7z.cloudfront.net/
unzip AGeNT.zip
rm AGeNT.zip
mkdir /AGeNT
mv LocatIt_v4.0.1.jar /AGeNT/LocatIt.jar
mv SurecallTrimmer_v4.0.1.jar /AGeNT/SurecallTrimmer.jar

# Install Pisces suite
cd /build
git clone https://github.com/Illumina/Pisces.git Pisces
cd Pisces
mkdir /Pisces/
cp -r binaries/5.2.7.47 /Pisces/
cd /Pisces/5.2.7.47
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done

# Clean up install
cd /
apt-get remove -y \
  build-essential \
  ca-certificates \
  gcc \
  git \
  libbz2-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  unzip \
  wget \
  zlib1g-dev
apt-get clean
rm -rf /build /var/lib/apt/lists/* /opt/get-pip.py