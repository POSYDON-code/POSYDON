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

# Version of POSYDON to install (release tag without the leading "v"),
# passed at build time by the CI workflow. Pinning guarantees the image
# contains the exact release version even if a newer one exists on the channel.
ARG POSYDON_VERSION

# Create a dedicated conda environment for POSYDON
RUN conda create --name posydon_env python=3.11 -y && \
    conda clean --all -f -y

# Install POSYDON from the posydon conda channel
RUN conda run -n posydon_env \
        conda install -c posydon -c conda-forge \
            posydon==${POSYDON_VERSION} \
            ipykernel \
            -y && \
    conda clean --all -f -y

# Fail the build if the installed version does not match the requested release
RUN conda run -n posydon_env \
        python -c "import posydon; assert posydon.__version__ == '${POSYDON_VERSION}', f'installed posydon {posydon.__version__} != requested {POSYDON_VERSION}'"

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
