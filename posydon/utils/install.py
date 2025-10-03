import os
import subprocess
from pathlib import Path

HOME = Path.home()
MINICONDA_DIR = HOME / "miniconda"
MINICONDA_SCRIPT = HOME / "miniconda.sh"
ENV_NAME = "posydon_env"
DATA_DIR = "/projects/e33022/POSYDON-shared/data"
ENV_VARS_DIR = HOME / "env_vars"
DOTENV_FILE = ENV_VARS_DIR / f"{ENV_NAME}.env"
IPYCONFIG_PATH = HOME / ".ipython" / "profile_default" / "ipython_kernel_config.py"

def run(cmd, check=True, shell=True, **kwargs):
    """Run a shell command."""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=shell, check=check, **kwargs)

def install_miniconda():
    if not MINICONDA_DIR.exists():
        print("📦 Installing Miniconda...")
        run(f"curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output {MINICONDA_SCRIPT}")
        run(f"bash {MINICONDA_SCRIPT} -b -p {MINICONDA_DIR}")
        MINICONDA_SCRIPT.unlink()

def activate_conda():
    print("⚡ Activating all shells for conda usage...")
    os.environ["PATH"] = f"{MINICONDA_DIR / 'bin'}:{os.environ['PATH']}"
    run("conda init --all")
    bashrc = HOME / ".bashrc"
    if bashrc.exists():
        run(f"source {bashrc}", shell=True, executable="/bin/bash")

def create_env_and_kernel():
    print("🐍 Creating conda environment and IPython kernel with POSYDON installed...")
    run("conda --version")
    run("conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main")
    run("conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r")
    run("conda config --add channels posydon")
    
    run(f"conda create --name {ENV_NAME} -c conda-forge -c posydon python=3.11 git posydon ipykernel --yes")
    
    activate_cmd = f"conda activate {ENV_NAME} && python -c \"import posydon; print(posydon.__version__)\""
    result = subprocess.run(activate_cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
    version = result.stdout.strip()
    
    run(f"conda activate {ENV_NAME} && python -m ipykernel install --user --name {ENV_NAME} --display-name \"POSYDON ({version})\"", executable="/bin/bash")
    
    return version

def set_conda_env_vars(posydon_dir):
    run(f"conda env config vars set PATH_TO_POSYDON={posydon_dir} --name {ENV_NAME}")
    run(f"conda env config vars set PATH_TO_POSYDON_DATA={DATA_DIR} --name {ENV_NAME}")

def setup_dotenv(posydon_dir):
    ENV_VARS_DIR.mkdir(parents=True, exist_ok=True)
    if DOTENV_FILE.exists():
        DOTENV_FILE.unlink()
    with open(DOTENV_FILE, "w") as f:
        f.write(f"PATH_TO_POSYDON={posydon_dir}\n")
        f.write(f"PATH_TO_POSYDON_DATA={DATA_DIR}\n")

def update_bashrc():
    bashrc = HOME / ".bashrc"
    conda_cmd = f"conda activate {ENV_NAME}"
    if bashrc.exists():
        with open(bashrc, "r") as f:
            lines = f.read().splitlines()
        if conda_cmd not in lines:
            with open(bashrc, "a") as f:
                f.write(f"\n{conda_cmd}\n")

def setup_ipython_config(posydon_dir):
    if IPYCONFIG_PATH.exists():
        IPYCONFIG_PATH.unlink()
    run("ipython profile create")
    tex_path = "/software/texlive/2020/bin/x86_64-linux/"
    ipy_config_line = (
        'c.InteractiveShellApp.exec_lines = ['
        f'"import os; os.environ[\'PATH\'] += \':{tex_path}\'; '
        f'os.environ[\'PATH_TO_POSYDON\'] = \'{posydon_dir}\'; '
        f'os.environ[\'PATH_TO_POSYDON_DATA\'] = \'{DATA_DIR}\'"]'
    )
    with open(IPYCONFIG_PATH, "a") as f:
        f.write(ipy_config_line + "\n")

def get_posydon_dir():
    cmd = f"conda activate {ENV_NAME} && python -c \"import posydon; print(posydon.__file__.split('/posydon/')[0])\""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
    return result.stdout.strip()

def _run_install():
    install_miniconda()
    activate_conda()
    version = create_env_and_kernel()
    posydon_dir = get_posydon_dir()
    set_conda_env_vars(posydon_dir)
    setup_dotenv(posydon_dir)
    update_bashrc()
    setup_ipython_config(posydon_dir)
    
    print(f"\n✅ Successfully installed POSYDON version {version}")
