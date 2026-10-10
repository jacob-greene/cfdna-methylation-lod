# wgbstools v0.3.0 (github.com/nloyfer/wgbs_tools) with samtools and tabix.
# Build: docker build -t cfdna-wgbstools:0.3.0 pipeline/containers/wgbstools
FROM python:3.11-slim-bookworm
RUN apt-get update \
 && apt-get install -y --no-install-recommends g++ make curl ca-certificates samtools tabix procps \
 && rm -rf /var/lib/apt/lists/*
COPY --from=ghcr.io/astral-sh/uv:0.4.30 /uv /usr/local/bin/uv
RUN uv pip install --system --no-cache 'numpy>=1.24,<2' 'pandas>=2,<3' 'scipy>=1.10,<2'
ARG WGBS_SHA256=e43e7ddbf535152ba5ed82bcc7e810ebf4bea90b4e3834f45703afb0a15773c5
RUN curl -fsSL -o /tmp/wgbs.tar.gz https://github.com/nloyfer/wgbs_tools/archive/refs/tags/v0.3.0.tar.gz \
 && echo "${WGBS_SHA256}  /tmp/wgbs.tar.gz" | sha256sum -c - \
 && mkdir -p /opt/wgbstools && tar xzf /tmp/wgbs.tar.gz -C /opt/wgbstools --strip-components=1 \
 && cd /opt/wgbstools && python setup.py && rm /tmp/wgbs.tar.gz \
 && chmod -R a+rX /opt/wgbstools \
 && ln -s /opt/wgbstools/src/python/wgbs_tools.py /usr/local/bin/wgbstools
