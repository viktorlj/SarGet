FROM davidcaste/alpine-java-unlimited-jce:latest

LABEL \
  author="Viktor Ljungström" \
  description="Image with tools used for NF targeted seq analysis, based on Max Ulysses Sarek but based on Java enabled image" \
  maintainer="viktor.ljungstrom@igp.uu.se"

# Export PATH
ENV \
  PATH=/opt/bcftools/bin:/opt/htslib/bin:/opt/samtools/bin:$PATH

# Build container
ADD \
  build.sh /usr/bin/
RUN \
  build.sh

CMD ["/bin/bash"]