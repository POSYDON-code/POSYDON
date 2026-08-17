# ============================================================
# Custom Jupyter image based on jupyter/scipy-notebook
# Adds a dedicated POSYDON conda environment registered
# as a named Jupyter kernel ("python (POSYDON)").
# ============================================================
FROM quay.io/jupyter/scipy-notebook:latest

# Install MESA SDK dependencies as root:
#   binutils     - development tools (for SDKs prior to 23.7.2)
#   make         - dependency/compilation tool
#   perl         - scripting language
#   libx11-dev   - X11 windowing library + development headers
#   zlib1g-dev   - Z compression library + development headers
#   tcsh         - the C shell / derivatives
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        binutils \
        make \
        perl \
        libx11-dev \
        zlib1g-dev \
        tcsh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

# Create a dedicated conda environment for POSYDON
RUN conda create --name posydon_env python=3.11 -y && \
    conda clean --all -f -y

# Install POSYDON from the posydon conda channel.
# POSYDON_VERSION is pinned to the released version by the release workflow so
# the image always matches the release; retries tolerate the Anaconda.org
# indexing delay after the package is uploaded.
ARG POSYDON_VERSION=*
RUN for i in 1 2 3 4 5; do \
        if conda run -n posydon_env \
            conda install -c posydon -c conda-forge \
                "posydon=${POSYDON_VERSION}" ipykernel -y; then \
            conda clean --all -f -y; \
            exit 0; \
        fi; \
        sleep 30; \
    done; \
    exit 1

# Register the posydon_env as a Jupyter kernel so it appears
# in the JupyterLab launcher and kernel selector as "Python (POSYDON)".
RUN conda run -n posydon_env \
        python -m ipykernel install \
            --prefix /opt/conda \
            --name posydon \
            --display-name "Python (POSYDON)"

# Set PATH_TO_POSYDON to installed conda version
ENV PATH_TO_POSYDON=/opt/conda/envs/posydon_env/lib/python3.11/site-packages

# Can be loaded from outside the container to point to a mounted data volume.
ENV PATH_TO_POSYDON_DATA=/home/${NB_USER}/data/POSYDON_GRIDS

USER root
RUN mkdir -p /etc/jupyter
COPY jupyter_server_config.py /etc/jupyter/jupyter_server_config.py
RUN rm -rf /opt/conda/share/jupyter/kernels/python3 && \
    fix-permissions /opt/conda/envs/posydon_env
USER ${NB_UID}
