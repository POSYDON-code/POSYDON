# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - Development Branch

### Added

#### Major New Features
- **New Command-Line Tools**: Complete suite of new command-line utilities for enhanced workflow management:
  - `compress-mesa`: Utility for compressing MESA files and managing storage
  - `get-posydon-data`: Tool for downloading and managing POSYDON data files
  - `posydon-popsyn`: Population synthesis execution tool
  - `posydon-run-grid`: Grid execution and management utility
  - `posydon-run-pipeline`: Pipeline execution tool for streamlined workflows
  - `posydon-setup-grid`: Grid setup and configuration utility
  - `posydon-setup-pipeline`: Pipeline setup and configuration tool

#### Documentation Improvements
- **Complete Documentation Restructure**: Comprehensive reorganization with new `_source` directory structure
- **New Documentation Components**:
  - `components-overview`: Detailed overview of POSYDON components and architecture
  - `contact-support`: User support and contact information
  - `contributing`: Guidelines for contributing to the project
  - `getting-started`: Streamlined getting started guide
  - `introduction-acknowledgements`: Project introduction and acknowledgements
  - `troubleshooting-faqs`: Comprehensive troubleshooting guide and FAQ
  - `tutorials-examples`: Extended tutorials and example workflows
- **Enhanced API Documentation**: Updated API documentation structure with improved navigation
- **New Population Parameter Guides**: Detailed guides for population synthesis parameters

#### Infrastructure and CI/CD
- **GitHub Actions Workflows**: Complete CI/CD pipeline implementation
  - `continuous_integration.yml`: Automated testing and integration
  - `deploy-github-pages-development.yml`: Development documentation deployment
  - `deploy-github-pages-release.yml`: Release documentation deployment
  - `install_extras.yml`: Extended installation testing
  - `publish-to-anaconda.yml`: Automated conda package publishing
- **Issue Templates**: New GitHub issue templates for better bug reporting and feature requests
  - `bug_report.md`: Structured bug report template
  - `feature_request.md`: Feature request template

#### Package and Dependencies
- **New Dependencies**: Enhanced functionality with additional packages
  - `hurry.filesize >= 0.9, <= 0.9`: File size utilities for better file management
  - `python-dotenv >= 1.0.0, <= 1.0.1`: Environment variable management
- **Enhanced Post-Processing Pipeline**: Improved pipeline capabilities with better error handling
- **Population Synthesis Tools**: New and improved population synthesis utilities

### Changed

#### Breaking Changes
- **Python Version Requirement**: Updated from Python 3.7 to Python 3.11
  - Minimum required version: `>=3.11, <3.12`
  - This is a breaking change that requires Python environment updates

#### Dependency Updates
- **Core Dependencies**: Updated version ranges for better stability and compatibility
  - `numpy >= 1.24.2, < 2.0.0` (updated from previous version)
  - `scipy >= 1.10.1, <= 1.14.1` (updated version range)
  - `iminuit >= 2.21.3, <= 2.30.1` (updated version range)
  - `configparser >= 5.3.0, <= 7.1.0` (updated version range)
  - `astropy >= 5.2.2, <= 6.1.6` (updated version range)
  - `pandas >= 2.0.0, <= 2.2.3` (updated version range)
  - `scikit-learn == 1.2.2` (pinned version)
  - `matplotlib >= 3.9.0, <= 3.9.2` (updated version range)
  - `matplotlib-label-lines >= 0.5.2, <= 0.7.0` (updated version range)
  - `h5py >= 3.8.0, <= 3.12.1` (updated version range)
  - `psutil >= 5.9.4, <= 6.1.0` (updated version range)
  - `tqdm >= 4.65.0, <= 4.67.0` (updated version range)
  - `tables >= 3.8.0, <= 3.10.1` (updated version range)
  - `progressbar2 >= 4.2.0, <= 4.5.0` (updated version range)

#### Package Structure
- **Enhanced Module Organization**: Improved organization of core modules
  - `active_learning`: Machine learning capabilities
  - `binary_evol`: Binary evolution modeling
  - `grids`: Grid management and interpolation
  - `interpolation`: Advanced interpolation methods
  - `popsyn`: Population synthesis tools
  - `utils`: Utility functions and helpers
  - `visualization`: Visualization and plotting tools

#### Build System
- **Updated Build Requirements**: 
  - `setuptools >= 76.0.0` (updated from previous version)
  - Enhanced build configuration with improved dependency management

### Fixed

#### Improvements and Bug Fixes
- **Enhanced Error Handling**: Better error messages and exception handling throughout the codebase
- **Improved File Management**: Enhanced file handling and compression utilities
- **Better Warning Systems**: More informative warnings and status messages
- **Conda Package Management**: Improved conda package configuration and dependencies

### Infrastructure

#### Legal and Licensing
- **License Update**: Updated copyright year to 2024
- **Project Metadata**: Updated project URLs and contact information
  - Homepage: https://posydon.org
  - Documentation: https://posydon.org/POSYDON
  - Repository: https://github.com/POSYDON-code/POSYDON.git
  - Issues: https://github.com/POSYDON-code/POSYDON/issues
  - Changelog: https://github.com/POSYDON-code/POSYDON/releases

#### Development Environment
- **Conda Environment**: Updated conda package specification with Python 3.11
- **Development Tools**: Enhanced development workflow with improved tooling
- **Testing Infrastructure**: Improved testing capabilities and coverage

---

## [v1.0.5] - Previous Release

### Notes
This changelog documents changes from v1.0.5 to the current development branch. For historical changes prior to v1.0.5, please refer to the git history or previous release notes.

---

## Migration Guide

### Upgrading from v1.0.5

#### Python Version
- **Required Action**: Update your Python environment to Python 3.11
- **Impact**: This is a breaking change that requires environment updates
- **Commands**:
  ```bash
  # Using conda
  conda create -n posydon-env python=3.11
  conda activate posydon-env
  
  # Using pyenv
  pyenv install 3.11.0
  pyenv local 3.11.0
  ```

#### Dependencies
- **Required Action**: Update all dependencies to new version ranges
- **Impact**: Some dependency versions may have breaking changes
- **Commands**:
  ```bash
  # Reinstall POSYDON with updated dependencies
  pip install --upgrade posydon
  
  # Or using conda
  conda update posydon
  ```

#### New Command-Line Tools
- **New Features**: Several new command-line tools are available
- **Impact**: Enhanced workflow capabilities
- **Usage**: See documentation for detailed usage instructions for each new tool

---

*This changelog is maintained by the POSYDON Collaboration. For questions or issues, please visit our [GitHub Issues](https://github.com/POSYDON-code/POSYDON/issues) page.*