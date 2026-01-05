# =============================================================================
# Stranger Container
# =============================================================================
# STR Annotation Tool
# https://github.com/moonso/stranger
#
# TODO: Customize this Dockerfile according to your environment needs
# =============================================================================

FROM python:3.11-slim

# Install stranger
RUN pip install --no-cache-dir stranger==0.8.0

# Install bcftools for VCF operations
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    tabix \
&& rm -rf /var/lib/apt/lists/*

# Copy custom scripts or configurations
# COPY bin/ /usr/local/bin/

# Default command
# CMD ["stranger"]
