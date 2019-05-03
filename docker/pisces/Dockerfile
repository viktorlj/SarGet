FROM microsoft/dotnet:2.1-sdk

LABEL \
  author="Viktor Ljungström" \
  description="Image with tools used for NF targeted seq analysis, dotnet images for Pisces" \
  maintainer="viktor.ljungstrom@igp.uu.se"

# Build container
ADD \
  build.sh /usr/bin/
RUN \
  build.sh

CMD ["/bin/bash"]