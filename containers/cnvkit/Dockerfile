# =============================================================================
# CNVkit Container
# =============================================================================
# Copy Number Variation detection and analysis
# https://github.com/etal/cnvkit
#
# TODO: Customize this Dockerfile according to your environment needs
# =============================================================================

FROM python:3.11-slim

# Install CNVkit
RUN pip install --no-cache-dir cnvkit==0.9.10

# Install additional tools for CNV analysis
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    tabix \
    samtools \
&& rm -rf /var/lib/apt/lists/*

# Copy custom scripts or configurations
# COPY bin/ /usr/local/bin/

# Default command
# CMD ["cnvkit"]
