#!/bin/bash

################################################################################
# POSYDON Installation Script
#
# This script automates the installation of POSYDON from the conda channel.
# It performs system checks, collects user inputs, creates a conda environment,
# installs POSYDON, configures environment variables, and downloads data.
#
# Quick Install (one-liner):
#   bash <(curl -sSL https://raw.githubusercontent.com/POSYDON-code/POSYDON/main/install.sh)
#
# Usage: bash install.sh
#
# Requirements:
#   - Conda/Anaconda/Miniconda installed (version >= 23.1.0 recommended)
#   - Internet connection for package and data downloads
#   - Sufficient disk space for selected datasets
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Global variables to store user inputs
ENV_NAME=""
DATA_PATH=""
DATASET_CHOICE=""
PATH_TO_POSYDON=""

################################################################################
# Helper Functions
################################################################################

print_header() {
    echo ""
    echo -e "${CYAN}${BOLD}════════════════════════════════════════════════════════════════${NC}"
    echo -e "${CYAN}${BOLD}  $1${NC}"
    echo -e "${CYAN}${BOLD}════════════════════════════════════════════════════════════════${NC}"
    echo ""
}

print_step() {
    echo -e "${BLUE}${BOLD}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}"
}

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to get conda version
get_conda_version() {
    conda --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' || echo "0.0.0"
}

# Function to compare versions
version_ge() {
    printf '%s\n%s' "$2" "$1" | sort -V -C
}

################################################################################
# Phase 1: Validation - System Prerequisites Check
################################################################################

check_prerequisites() {
    print_header "Phase 1: System Prerequisites Check"

    # Check if conda is installed
    print_step "Checking for conda installation..."
    if ! command_exists conda; then
        print_error "Conda is not installed or not in PATH!"
        echo ""
        echo -e "${BOLD}Conda/Anaconda/Miniconda is required to install POSYDON.${NC}"
        echo ""
        echo "Please install one of the following:"
        echo ""
        echo "  • Miniconda (recommended, minimal installer):"
        echo "    https://docs.conda.io/en/latest/miniconda.html"
        echo ""
        echo "  • Anaconda (full distribution with many packages):"
        echo "    https://www.anaconda.com/products/distribution"
        echo ""
        echo "After installation, restart your terminal and run this script again."
        echo ""
        exit 1
    fi
    print_success "Conda found: $(command -v conda)"

    # Check conda version
    print_step "Checking conda version..."
    CONDA_VERSION=$(get_conda_version)
    print_info "Detected conda version: $CONDA_VERSION"

    if ! version_ge "$CONDA_VERSION" "23.1.0"; then
        print_warning "Your conda version ($CONDA_VERSION) is older than the recommended version (23.1.0)."
        print_warning "Installation may be very slow. Consider updating conda with:"
        echo "    conda update -n base conda"
        echo ""
        read -p "Continue anyway? (y/n): " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Installation cancelled."
            exit 0
        fi
    else
        print_success "Conda version is sufficient (>= 23.1.0)"
    fi

    # Check for libmamba solver (optional but recommended)
    print_step "Checking conda solver configuration..."
    SOLVER=$(conda config --show solver 2>/dev/null | grep -oE 'libmamba|classic' || echo "classic")
    if [[ "$SOLVER" == "libmamba" ]]; then
        print_success "Using libmamba solver (fast dependency resolution)"
    else
        print_warning "Using classic solver (slower dependency resolution)"
        print_info "For faster installation, consider configuring libmamba solver:"
        echo "    conda install -n base conda-libmamba-solver"
        echo "    conda config --set solver libmamba"
        echo ""
    fi

    echo ""
    print_success "All prerequisite checks passed!"
    echo ""
}

################################################################################
# Phase 2: Input Collection - Gather User Inputs
################################################################################

collect_environment_name() {
    print_step "Enter conda environment name"
    echo ""
    echo "This will create a new conda environment for POSYDON."
    echo -n "Environment name [default: posydon-env]: "
    read -r ENV_NAME

    # Use default if empty
    if [[ -z "$ENV_NAME" ]]; then
        ENV_NAME="posydon-env"
    fi

    # Validate environment name (alphanumeric, hyphens, underscores only)
    if [[ ! "$ENV_NAME" =~ ^[a-zA-Z0-9_-]+$ ]]; then
        print_error "Invalid environment name! Use only letters, numbers, hyphens, and underscores."
        echo ""
        exit 1
    fi

    # Check if environment already exists (more robust check)
    print_step "Checking if environment '${ENV_NAME}' already exists..."
    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        print_error "Environment '${ENV_NAME}' already exists!"
        echo ""
        echo "Options:"
        echo "  1. Re-run this script and choose a different name"
        echo "  2. Remove existing environment with: conda env remove -n ${ENV_NAME}"
        echo ""
        exit 1
    fi

    print_success "Environment name '${ENV_NAME}' is available"
    echo ""
}

collect_data_path() {
    print_step "Enter POSYDON data directory path"
    echo ""
    echo "This directory will store POSYDON datasets (can be very large, 10-100+ GB)."
    echo "The path will be set as PATH_TO_POSYDON_DATA environment variable."
    echo ""

    # Check if PATH_TO_POSYDON_DATA is already set in the environment
    if [[ -n "${PATH_TO_POSYDON_DATA:-}" ]]; then
        print_info "PATH_TO_POSYDON_DATA is already set to: ${PATH_TO_POSYDON_DATA}"
        echo -n "Use this path? (y/n): "
        read -r -n 1 USE_EXISTING
        echo ""
        if [[ $USE_EXISTING =~ ^[Yy]$ ]]; then
            DATA_PATH="$PATH_TO_POSYDON_DATA"
            print_success "Using existing path: ${DATA_PATH}"
            echo ""
            return 0
        fi
        echo ""
    fi

    while true; do
        echo -n "Data directory path: "
        read -r DATA_PATH

        # Expand tilde to home directory
        DATA_PATH="${DATA_PATH/#\~/$HOME}"

        if [[ -z "$DATA_PATH" ]]; then
            print_error "Path cannot be empty!"
            continue
        fi

        # Check if directory exists
        if [[ -d "$DATA_PATH" ]]; then
            print_success "Directory exists: ${DATA_PATH}"

            # Check write permissions
            if [[ -w "$DATA_PATH" ]]; then
                print_success "Write permission confirmed"
                break
            else
                print_error "No write permission for directory: ${DATA_PATH}"
                echo "Please choose a directory where you have write access."
                echo ""
            fi
        else
            # Ask to create directory
            echo ""
            print_warning "Directory does not exist: ${DATA_PATH}"
            echo -n "Create this directory? (y/n): "
            read -r -n 1 CREATE_DIR
            echo ""

            if [[ $CREATE_DIR =~ ^[Yy]$ ]]; then
                if mkdir -p "$DATA_PATH" 2>/dev/null; then
                    print_success "Directory created: ${DATA_PATH}"
                    break
                else
                    print_error "Failed to create directory. Check parent directory permissions."
                    echo ""
                fi
            else
                echo "Please enter a different path."
                echo ""
            fi
        fi
    done

    echo ""
}

collect_dataset_choice() {
    print_step "Select POSYDON dataset to download"
    echo ""
    echo "Available complete datasets (include auxiliary data):"
    echo ""
    echo -e "  ${BOLD}Complete Sets:${NC}"
    echo -e "    1) DR2          - Full DR2 dataset (~103 GB, all 8 metallicities + auxiliary) ${GREEN}[Recommended]${NC}"
    echo -e "    2) DR2_1Zsun    - Solar metallicity (~14 GB, 1 Zsun + auxiliary)"
    echo "    3) DR2_2Zsun    - Twice solar metallicity (~14 GB, 2 Zsun + auxiliary)"
    echo "    4) DR2_0.45Zsun - Sub-solar metallicity (~14 GB, 0.45 Zsun + auxiliary)"
    echo "    5) DR2_0.2Zsun  - Low metallicity (~14 GB, 0.2 Zsun + auxiliary)"
    echo "    6) DR2_0.1Zsun  - Low metallicity (~14 GB, 0.1 Zsun + auxiliary)"
    echo "    7) DR2_0.01Zsun - Very low metallicity (~14 GB, 0.01 Zsun + auxiliary)"
    echo "    8) DR2_1e-3Zsun - Extremely low metallicity (~14 GB, 10^-3 Zsun + auxiliary)"
    echo "    9) DR2_1e-4Zsun - Extremely low metallicity (~14 GB, 10^-4 Zsun + auxiliary)"
    echo ""
    echo -e "  ${BOLD}Legacy Sets:${NC}"
    echo "   10) DR1          - Data Release 1 (~40 GB, solar metallicity)"
    echo "   11) DR1.1        - DR1 + super-Eddington (~40 GB, solar metallicity)"
    echo ""
    echo -e "  ${BOLD}Minimal:${NC}"
    echo "   12) auxiliary    - Auxiliary data only (~4 GB, required for all runs)"
    echo ""
    echo "   13) Skip         - Skip data download (install manually later)"
    echo ""

    while true; do
        echo -n "Enter choice [1-13]: "
        read -r CHOICE

        case $CHOICE in
            1) DATASET_CHOICE="DR2"; DATASET_SIZE="~103 GB"; break ;;
            2) DATASET_CHOICE="DR2_1Zsun"; DATASET_SIZE="~14 GB"; break ;;
            3) DATASET_CHOICE="DR2_2Zsun"; DATASET_SIZE="~14 GB"; break ;;
            4) DATASET_CHOICE="DR2_0.45Zsun"; DATASET_SIZE="~14 GB"; break ;;
            5) DATASET_CHOICE="DR2_0.2Zsun"; DATASET_SIZE="~14 GB"; break ;;
            6) DATASET_CHOICE="DR2_0.1Zsun"; DATASET_SIZE="~14 GB"; break ;;
            7) DATASET_CHOICE="DR2_0.01Zsun"; DATASET_SIZE="~14 GB"; break ;;
            8) DATASET_CHOICE="DR2_1e-3Zsun"; DATASET_SIZE="~14 GB"; break ;;
            9) DATASET_CHOICE="DR2_1e-4Zsun"; DATASET_SIZE="~14 GB"; break ;;
            10) DATASET_CHOICE="DR1"; DATASET_SIZE="~40 GB"; break ;;
            11) DATASET_CHOICE="DR1.1"; DATASET_SIZE="~40 GB"; break ;;
            12) DATASET_CHOICE="auxiliary"; DATASET_SIZE="~4 GB"; break ;;
            13) DATASET_CHOICE="skip"; DATASET_SIZE="N/A"; break ;;
            *) print_error "Invalid choice. Please enter a number between 1 and 13."; continue ;;
        esac
    done

    if [[ "$DATASET_CHOICE" != "skip" ]]; then
        print_success "Selected dataset: ${DATASET_CHOICE} (${DATASET_SIZE})"
        print_warning "Ensure you have at least ${DATASET_SIZE} of free disk space in ${DATA_PATH}"
    else
        print_info "Data download will be skipped"
    fi
    echo ""
}

confirm_inputs() {
    print_header "Phase 2: Configuration Summary"
    echo ""
    echo "Please review your configuration:"
    echo ""
    echo -e "  ${BOLD}Conda Environment:${NC}     ${ENV_NAME}"
    echo -e "  ${BOLD}Data Directory:${NC}        ${DATA_PATH}"
    echo -e "  ${BOLD}Dataset:${NC}               ${DATASET_CHOICE}"
    if [[ "$DATASET_CHOICE" != "skip" ]]; then
        echo -e "  ${BOLD}Estimated Size:${NC}        ${DATASET_SIZE}"
    fi
    echo ""

    while true; do
        echo -n "Proceed with installation? (y/n): "
        read -r -n 1 CONFIRM
        echo ""

        if [[ $CONFIRM =~ ^[Yy]$ ]]; then
            print_success "Configuration confirmed!"
            echo ""
            return 0
        elif [[ $CONFIRM =~ ^[Nn]$ ]]; then
            echo "Installation cancelled by user."
            exit 0
        else
            print_error "Please enter 'y' or 'n'"
        fi
    done
}

collect_user_inputs() {
    print_header "Phase 2: User Input Collection"
    echo ""

    collect_environment_name
    collect_data_path
    collect_dataset_choice
    confirm_inputs
}

################################################################################
# Phase 3: Environment Creation
################################################################################

create_conda_environment() {
    print_header "Phase 3: Conda Environment Creation"

    print_step "Creating conda environment '${ENV_NAME}' with Python 3.11..."
    echo ""

    if conda create --name "$ENV_NAME" python=3.11 -y; then
        print_success "Conda environment created successfully"
    else
        print_error "Failed to create conda environment"
        exit 1
    fi
    echo ""

    print_step "Adding conda-forge channel..."
    if conda config --add channels conda-forge 2>/dev/null; then
        print_success "conda-forge channel configured"
    else
        print_warning "conda-forge channel may already be configured"
    fi
    echo ""

    print_success "Environment setup complete"
    echo ""
}

################################################################################
# Phase 4: POSYDON Installation
################################################################################

install_posydon() {
    print_header "Phase 4: POSYDON Installation"

    print_step "Installing POSYDON from conda channel..."
    print_info "This may take several minutes to resolve dependencies..."
    echo ""

    # Source conda to ensure conda activate works in script
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"

    # Activate environment
    if conda activate "$ENV_NAME"; then
        print_success "Environment '${ENV_NAME}' activated"
    else
        print_error "Failed to activate environment"
        exit 1
    fi
    echo ""

    # Install POSYDON
    print_step "Running: conda install -c posydon posydon -y"
    echo ""
    if conda install -c posydon posydon -y; then
        print_success "POSYDON installed successfully"
    else
        print_error "Failed to install POSYDON"
        exit 1
    fi
    echo ""

    # Verify installation and discover PATH_TO_POSYDON
    print_step "Verifying POSYDON installation..."
    PATH_TO_POSYDON=$(python -c "import posydon, os; print(os.path.dirname(os.path.dirname(posydon.__file__)))" 2>/dev/null)

    if [[ -n "$PATH_TO_POSYDON" ]] && [[ -d "$PATH_TO_POSYDON" ]]; then
        print_success "POSYDON installed at: ${PATH_TO_POSYDON}"
    else
        print_error "Failed to verify POSYDON installation"
        exit 1
    fi
    echo ""

    print_success "POSYDON installation complete"
    echo ""
}

################################################################################
# Phase 5: Environment Variable Configuration
################################################################################

configure_environment_variables() {
    print_header "Phase 5: Environment Variable Configuration"

    print_step "Setting PATH_TO_POSYDON environment variable..."
    if conda env config vars set PATH_TO_POSYDON="$PATH_TO_POSYDON" -n "$ENV_NAME"; then
        print_success "PATH_TO_POSYDON = ${PATH_TO_POSYDON}"
    else
        print_error "Failed to set PATH_TO_POSYDON"
        exit 1
    fi
    echo ""

    print_step "Setting PATH_TO_POSYDON_DATA environment variable..."
    if conda env config vars set PATH_TO_POSYDON_DATA="$DATA_PATH" -n "$ENV_NAME"; then
        print_success "PATH_TO_POSYDON_DATA = ${DATA_PATH}"
    else
        print_error "Failed to set PATH_TO_POSYDON_DATA"
        exit 1
    fi
    echo ""

    print_info "Environment variables will be available after reactivating the environment"
    echo ""

    print_success "Environment variables configured"
    echo ""
}

################################################################################
# Phase 6: Data Download
################################################################################

download_data() {
    print_header "Phase 6: Data Download"

    if [[ "$DATASET_CHOICE" == "skip" ]]; then
        print_info "Data download skipped as requested"
        echo ""
        return 0
    fi

    print_step "Reactivating environment to load environment variables..."
    conda deactivate
    conda activate "$ENV_NAME"
    print_success "Environment reactivated"
    echo ""

    print_step "Downloading POSYDON dataset: ${DATASET_CHOICE}"
    print_info "This may take a considerable amount of time depending on your connection..."
    print_info "Download size: ${DATASET_SIZE}"
    echo ""

    if get-posydon-data "$DATASET_CHOICE"; then
        print_success "Dataset downloaded successfully!"
    else
        print_error "Failed to download dataset"
        print_warning "You can manually download the data later with:"
        echo "    conda activate ${ENV_NAME}"
        echo "    get-posydon-data ${DATASET_CHOICE}"
        echo ""
        echo "Or list available datasets with:"
        echo "    get-posydon-data -l complete"
        exit 1
    fi
    echo ""

    print_success "Data download complete"
    echo ""
}

################################################################################
# Phase 7: Final Instructions
################################################################################

show_final_instructions() {
    print_header "Phase 7: Installation Complete! 🎉"

    print_success "POSYDON has been successfully installed!"
    echo ""
    echo -e "${BOLD}Installation Summary:${NC}"
    echo -e "  • Environment name:      ${GREEN}${ENV_NAME}${NC}"
    echo -e "  • Python version:        ${GREEN}3.11${NC}"
    echo -e "  • POSYDON location:      ${GREEN}${PATH_TO_POSYDON}${NC}"
    echo -e "  • Data location:         ${GREEN}${DATA_PATH}${NC}"
    if [[ "$DATASET_CHOICE" != "skip" ]]; then
        echo -e "  • Dataset installed:     ${GREEN}${DATASET_CHOICE}${NC}"
    else
        echo -e "  • Dataset installed:     ${YELLOW}None (skipped)${NC}"
    fi
    echo ""

    echo -e "${BOLD}Next Steps:${NC}"
    echo ""
    echo -e "1. ${BOLD}Activate the POSYDON environment:${NC}"
    echo -e "   ${CYAN}conda activate ${ENV_NAME}${NC}"
    echo ""

    if [[ "$DATASET_CHOICE" == "skip" ]]; then
        echo -e "2. ${BOLD}Download POSYDON data:${NC}"
        echo "   First, list available datasets:"
        echo -e "   ${CYAN}get-posydon-data -l complete${NC}"
        echo ""
        echo "   Then download your chosen dataset (e.g., for solar metallicity):"
        echo -e "   ${CYAN}get-posydon-data DR2_1Zsun${NC}"
        echo ""
        echo -e "3. ${BOLD}Start using POSYDON:${NC}"
    else
        echo -e "2. ${BOLD}Start using POSYDON:${NC}"
    fi

    echo -e "${BOLD}Environment Variables:${NC}"
    echo "  The following environment variables are set in your conda environment:"
    echo -e "  • ${CYAN}PATH_TO_POSYDON${NC} = ${PATH_TO_POSYDON}"
    echo -e "  • ${CYAN}PATH_TO_POSYDON_DATA${NC} = ${DATA_PATH}"
    echo ""
    echo "  These are automatically loaded when you activate the environment."
    echo ""

    echo -e "${BOLD}Documentation:${NC}"
    echo -e "  • POSYDON Documentation: ${CYAN}https://posydon.org/POSYDON/latest/${NC}"
    echo ""

    print_success "Installation complete! Happy POSYDON-ing! 🌟"
    echo ""
}

################################################################################
# Main Execution
################################################################################

main() {
    clear

    echo ""
    echo -e "${BOLD}${CYAN}"
    echo "╔═══════════════════════════════════════════════════════════╗"
    echo "║                                                           ║"
    echo "║                POSYDON Installation Script                ║"
    echo "║                                                           ║"
    echo "║    Population Synthesis with Detailed Binary Evolution    ║"
    echo "║                                                           ║"
    echo "╚═══════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo ""
    echo "This script will guide you through installing POSYDON and setting up"
    echo "your environment."
    echo ""

    # Execute installation phases
    check_prerequisites
    collect_user_inputs
    create_conda_environment
    install_posydon
    configure_environment_variables
    download_data
    show_final_instructions
}

# Run main function
main
