ARG BASE_IMAGE=kga-splice:v6.3
FROM ${BASE_IMAGE}

# Keep regtools isolated from the existing analysis and annotation environments.
RUN conda create -y -n clinv1-regtools --override-channels \
        -c conda-forge -c bioconda "regtools=1.0.0" && \
    conda clean -y --tarballs

ENV PATH="/opt/conda/envs/clinv1-regtools/bin:/opt/conda/envs/KGA_Splice/bin:/opt/conda/bin:${PATH}" \
    VEP_CACHE_DIR="/refs/vep" \
    VEP_ASSEMBLY="GRCh38"

WORKDIR /opt/clinv1
COPY scripts/ /opt/clinv1/scripts/
COPY tests/ /opt/clinv1/tests/

# Fail the image build if the inherited environments or bundled workflow are unusable.
RUN command -v Rscript && command -v samtools && command -v sambamba && command -v bedtools && command -v gatk && \
    command -v regtools && command -v conda && \
    Rscript -e 'stopifnot(requireNamespace("xml2", quietly=TRUE))' && \
    conda run --no-capture-output -n vep_env vep --help > /dev/null && \
    conda run --no-capture-output -n spliceai_env spliceai -h > /dev/null && \
    conda run --no-capture-output -n KGA_Splice sashimipy --help > /dev/null && \
    Rscript tests/test_core.R && \
    Rscript tests/test_report.R && \
    Rscript tests/test_sai.R && \
    Rscript tests/test_report_pipeline.R && \
    Rscript tests/test_reference_pool.R && \
    Rscript tests/test_reference_bank.R && \
    bash tests/test_mock_pipeline.sh && \
    bash tests/test_stranded_coverage.sh && \
    bash tests/test_region_extraction.sh && \
    bash tests/test_reference_workflows.sh

CMD ["bash", "--noprofile", "--norc"]
