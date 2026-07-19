# ============================================================
# Custom Jupyter image based on jupyter/scipy-notebook
# Adds a dedicated POSYDON conda environment registered
# as a named Jupyter kernel ("python (POSYDON)").
# ============================================================
FROM quay.io/jupyter/scipy-notebook:latest

USER ${NB_UID}

# Create a dedicated conda environment for POSYDON
RUN conda create --name posydon_env python=3.11 -y && \
    conda clean --all -f -y

# Install POSYDON from the posydon conda channel
RUN conda run -n posydon_env \
        conda install -c posydon -c conda-forge \
            posydon \
            ipykernel \
            -y && \
    conda clean --all -f -y

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
