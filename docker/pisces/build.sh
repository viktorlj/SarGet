mkdir /build

# Install Pisces suite
cd /build
git clone https://github.com/Illumina/Pisces.git Pisces
cd Pisces
mkdir /Pisces/
cp -r binaries/5.2.10.49 /Pisces/
cd /Pisces/5.2.10.49
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done