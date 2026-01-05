# Container Registry

This directory contains Dockerfiles for all tools used in the pipeline.

## Status Overview

| Container | modules.config Image | Dockerfile | Status |
|-----------|---------------------|------------|--------|
| Mapping | `ghcr.io/pzweuj/mapping:2025dec` | mapping/Dockerfile | ✓ Ready |
| GATK | `broadinstitute/gatk:4.6.2.0` | gatk/Dockerfile | ✓ Ready |
| DeepVariant | `google/deepvariant:1.6.0` | deepvariant/Dockerfile | ✓ Ready |
| **GLNexus** | `ghcr.io/your-org/glnexus:v1.4.1` | glnexus/Dockerfile | ⚠️ Placeholder |
| VEP | `ensemblorg/ensembl-vep:release_112.0` | vep/Dockerfile | ✓ Ready |
| ExpansionHunter | `clinicalgenomics/expansionhunter:5.0.0` | expansionhunter/Dockerfile | ✓ Ready |
| **Stranger** | `ghcr.io/your-org/stranger:0.8.0` | stranger/Dockerfile | ⚠️ Placeholder |
| **CNVkit** | `ghcr.io/your-org/cnvkit:0.9.10` | cnvkit/Dockerfile | ⚠️ Placeholder |
| GenMod | `ghcr.io/pzweuj/mapping:2025dec` | genmod/Dockerfile | ✓ Ready |
| WhatsHap | `ghcr.io/pzweuj/mapping:2025dec` | whatshap/Dockerfile | ✓ Ready |
| PLINK2 | `ghcr.io/pzweuj/mapping:2025dec` | plink2/Dockerfile | ✓ Ready |

## Pending Actions (TODO)

The following containers need to be built and pushed to your registry:

```bash
# Build GLNexus
cd containers/glnexus
docker build -t ghcr.io/your-org/glnexus:v1.4.1 .
docker push ghcr.io/your-org/glnexus:v1.4.1

# Build Stranger
cd containers/stranger
docker build -t ghcr.io/your-org/stranger:0.8.0 .
docker push ghcr.io/your-org/stranger:0.8.0

# Build CNVkit
cd containers/cnvkit
docker build -t ghcr.io/your-org/cnvkit:0.9.10 .
docker push ghcr.io/your-org/cnvkit:0.9.10
```

Then update `conf/modules.config` to use your registry images.

## Building Containers

```bash
# Build a specific container
docker build -t ghcr.io/your-org/glnexus:v1.4.1 ./glnexus

# Build all containers
for dir in */; do
    name=$(basename "$dir")
    docker build -t ghcr.io/your-org/${name}:latest ./${name}
done
```

## Pushing to Registry

```bash
# Login to registry
docker login ghcr.io

# Push a specific container
docker push ghcr.io/your-org/glnexus:v1.4.1
```

## Notes

- Placeholder Dockerfiles (marked with ⚠️) need to be customized according to your HPC environment
- Some containers use official images from upstream projects as base images
- For HPC environments, consider building Singularity/Apptainer images from Docker images
