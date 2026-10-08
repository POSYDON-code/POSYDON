
import copy
import os
import re
from abc import ABC, abstractmethod


class InlistSection:
    """Represents a single namelist section (e.g., &star_job, &controls)"""
    def __init__(self, name, parameters=None):
        self.name = name
        self.parameters = parameters or {}

    def __repr__(self):
        """Return a string representation of the section"""
        param_count = len(self.parameters)
        if param_count == 0:
            return f"InlistSection(name='{self.name}', parameters=0)"

        # Show first few parameters as preview
        preview_keys = list(self.parameters.keys())[:3]
        preview = ", ".join(preview_keys)
        if param_count > 3:
            preview += ", ..."

        return f"InlistSection(name='{self.name}', parameters={param_count}: [{preview}])"

    def merge(self, other):
        """Merge another section into this one (later values override)"""
        if self.name != other.name:
            raise ValueError("Cannot merge sections with different names")
        self.parameters.update(other.parameters)
        return self

    def to_string(self):
        """Convert the section to a string representation"""
        lines = [f"&{self.name}"]
        for key, value in self.parameters.items():
            lines.append(f"    {key} = {value}")
        lines.append("/\n")
        return "\n".join(lines)

    def to_fortran(self) -> str:
        """Convert to Fortran namelist format"""
        lines = [f"&{self.name}"]
        for key, value in self.parameters.items():
            lines.append(f"\t{key} = {self._format_value(value)}")
        lines.append(f"/ ! end of {self.name} namelist\n")
        return "\n".join(lines)

    @staticmethod
    def _format_value(value) -> str:
        """Return value as-is (no formatting needed)"""
        if isinstance(value, bool):
            return '.true.' if value else '.false.'
        return value

    @staticmethod
    def _parse_value(value_str):
        """Return parsed value"""

        return value_str.strip()

    @classmethod
    def from_string(cls, name: str, content: str):
        """Parse a namelist section from string content"""
        parameters = {}

        # Split into lines and process
        lines = content.split('\n')
        for line in lines:
            # Remove comments
            if '!' in line:
                line = line[:line.index('!')]
            line = line.strip()

            # Skip empty lines, section headers, and end markers
            if not line or line.startswith('&') or line == '/':
                continue

            # Parse key = value
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                parameters[key] = cls._parse_value(value)

        return cls(name, parameters)

class Inlist:
    """Represents a complete inlist file with multiple sections"""
    def __init__(self, name, sections=None):
        object.__setattr__(self, 'name', name)
        object.__setattr__(self, 'sections', sections if sections is not None else {})

    def __repr__(self):
        """Return a string representation of the inlist"""
        section_count = len(self.sections)
        if section_count == 0:
            return f"Inlist(name='{self.name}', sections=0)"

        section_names = list(self.sections.keys())
        sections_str = ", ".join(section_names)

        return f"Inlist(name='{self.name}', sections={section_count}: [{sections_str}])"

    def __getattr__(self, name):
        """Allow attribute-style access to sections (e.g., inlist.controls)"""
        # Avoid infinite recursion by checking __dict__ first
        if name in ['name', 'sections']:
            return object.__getattribute__(self, name)

        # Try to find section by name
        sections = object.__getattribute__(self, 'sections')
        if name in sections:
            return sections[name]

        # If not found, raise AttributeError
        raise AttributeError(f"Inlist has no section '{name}'")

    def __setattr__(self, name, value):
        """Allow setting sections as attributes"""
        # Handle the standard attributes
        if name in ['name', 'sections']:
            object.__setattr__(self, name, value)
        # If value is an InlistSection, add it to sections
        elif isinstance(value, InlistSection):
            self.sections[name] = value
        else:
            # For other attributes, use normal behavior
            object.__setattr__(self, name, value)

    def add_section(self, section):
        """Add or merge a section into the inlist"""
        if section.name in self.sections:
            self.sections[section.name].merge(section)
        else:
            self.sections[section.name] = section

    def merge(self, other):
        """Merge another inlist into this one (later values override)"""
        merged = Inlist(self.name, copy.deepcopy(self.sections))
        for section in other.sections.values():
            merged.add_section(section)
        return merged

    def to_file(self, filepath):
        """Generate the complete inlist file content"""
        lines = []

        # Add sections
        # sort sections controls before star_job and binary_control before binary_job
        section_order = ['controls', 'star_job', 'binary_controls', 'binary_job']
        for section_name in section_order:
            if section_name in self.sections:
                lines.append(self.sections[section_name].to_fortran())

        for section_name, section in self.sections.items():
            if section_name not in section_order:
                lines.append(section.to_fortran())

        with open(filepath, 'w') as f:
            f.write("\n".join(lines))


    @classmethod
    def from_file(cls, filepath: str, name: str = None, section: str = None):
        """Read and parse an inlist file

        Args:
            filepath: Path to the inlist file
            name: Optional name for the inlist (defaults to filename)
            section: Optional section name if file has no &section markers (e.g., .defaults files)

        Returns:
            Inlist object with parsed sections
        """
        import os

        if name is None:
            name = os.path.basename(filepath)

        with open(filepath, 'r') as f:
            content = f.read()

        return cls.from_string(content, name, section=section)

    @classmethod
    def from_string(cls, content: str, name: str = "inlist", section: str = None):
        """Parse an inlist from string content

        This follows the same logic as clean_inlist_file for consistency.

        Args:
            content: String content of the inlist file
            name: Name for the inlist
            section: Optional section name if file has no &section markers (e.g., .defaults files)

        Returns:
            Inlist object with parsed sections
        """
        inlist = cls(name)

        # Clean inlist into nice list (matching clean_inlist_file logic)
        param_value_list = []
        for line in content.split('\n'):
            # Strip away all the comments and whitespace
            param_and_value = line.strip('\n').strip().split('!')[0].strip()
            # Check that this line actually has a parameter and value pair
            if param_and_value and ('=' in param_and_value or '&' in param_and_value):
                param_value_list.append(param_and_value)

        # Does this inlist have multiple sections?
        sections_found = {k: {} for k in param_value_list if '&' in k}

        if not sections_found:
            # MESA default files do not have sections,
            # because controls and jobs are in separate files
            if section:
                parameters = {}
                for item in param_value_list:
                    if '=' in item:
                        key, value = item.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        # Skip empty values like '', '.'
                        if value not in ["''", "'.'"]:
                            parameters[key] = InlistSection._parse_value(value)

                inlist.add_section(InlistSection(section, parameters))
        else:
            # Inlist has both job and controls sections marked with &
            current_section = None
            for item in param_value_list:
                if '&' in item:
                    # New section header
                    section_name = item.replace('&', '').strip()
                    current_section = section_name
                    sections_found[item] = {}
                elif '=' in item and current_section:
                    key, value = item.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    # Skip empty values like '', '.'
                    if value not in ["''", "'.'"]:
                        sections_found['&' + current_section][key] = InlistSection._parse_value(value)

            # Convert to InlistSection objects
            for section_key, params in sections_found.items():
                section_name = section_key.replace('&', '').strip()
                inlist.add_section(InlistSection(section_name, params))

        return inlist


class MESAInlists():
    """Handles MESA inlists for single and binary star evolution."""
    def __init__(self, path):
        self.path = path
        # .defaults files don't have &section markers, so we specify the section name
        controls_inlist = Inlist.from_file(f'{self.path}/star/controls.defaults', section='controls')
        star_job_inlist = Inlist.from_file(f'{self.path}/star/star_job.defaults', section='star_job')
        self.base_star_inlist = controls_inlist.merge(star_job_inlist)

        binary_job_inlist = Inlist.from_file(f'{self.path}/binary/binary_job.defaults', section='binary_job')
        binary_controls_inlist = Inlist.from_file(f'{self.path}/binary/binary_controls.defaults', section='binary_controls')
        self.base_binary_inlist = binary_job_inlist.merge(binary_controls_inlist)

        # Clean up default parameters (remove read_extra/inlist refs, replace num_x_ctrls placeholders)
        sections_to_clean = [
            self.base_star_inlist.controls,
            self.base_star_inlist.star_job,
            self.base_binary_inlist.binary_controls,
            self.base_binary_inlist.binary_job,
        ]
        for section in sections_to_clean:
            section.parameters = self._clean_parameters(section.parameters)

    def _clean_parameters(self, params_dict):
        """Clean parameters by removing read_extra/inlist references and replacing num_x_ctrls.

        Parameters
        ----------
            params_dict: dict
                Dictionary of parameters to clean

        Returns
        -------
            dict
                Cleaned parameters dictionary
        """
        # Remove read_extra and inlist references
        cleaned = {k: v for k, v in params_dict.items()
                   if not any(substring in k for substring in ['read_extra', 'inlist'])}

        # Replace num_x_ctrls with actual index (placeholder in MESA defaults)
        keys_to_replace = {k: k.replace('num_x_ctrls', '1')
                           for k in cleaned.keys() if 'num_x_ctrls' in k}
        for old_key, new_key in keys_to_replace.items():
            cleaned[new_key] = cleaned.pop(old_key)

        return cleaned

    def __repr__(self):
        return (f"MESAInlists(path='{self.path}', "
                f"base_star_inlist={self.base_star_inlist}, "
                f"base_binary_inlist={self.base_binary_inlist})")

class InlistManager:
    """Manages multiple inlists for different evolution steps/phases in MESA. Each entry in the lists corresponds to a different step/phase in the evolution sequence."""

    def __init__(self,):
        self.binary_inlists = []
        self.binary_star1_inlists = []
        self.binary_star2_inlists = []
        self.star1_inlists = []
        self.star2_inlists = []

    def keys(self):
        return ['binary_inlists', 'binary_star1_inlists', 'binary_star2_inlists', 'star1_inlists', 'star2_inlists']

    def __getitem__(self, key):
        if key in self.keys():
            return getattr(self, key)
        else:
            raise KeyError(f"InlistManager has no key '{key}'")


    def append_binary_inlist(self, inlist):
        self.binary_inlists.append(inlist)

    def append_binary_star1_inlist(self, inlist):
        self.binary_star1_inlists.append(inlist)

    def append_binary_star2_inlist(self, inlist):
        self.binary_star2_inlists.append(inlist)

    def append_star1_inlist(self, inlist):
        self.star1_inlists.append(inlist)

    def append_star2_inlist(self, inlist):
        self.star2_inlists.append(inlist)

    def __repr__(self):
        return (f"InlistManager(\nbinary_inlists={len(self.binary_inlists)}, \n"
                f"binary_star1_inlists={len(self.binary_star1_inlists)}, \n"
                f"binary_star2_inlists={len(self.binary_star2_inlists)}, \n"
                f"star1_inlists={len(self.star1_inlists)}, \n"
                f"star2_inlists={len(self.star2_inlists)})")

    def _write_inlist_group(self, inlists, output_path, single_name, step_name_pattern):
        """Helper to write a group of inlists with consistent naming logic.

        Parameters
        ----------
        inlists: List of Inlist objects
            List of inlists to write
        output_path: str
            Directory path to write to
        single_name: str
            Filename if only one inlist (e.g., 'inlist_project')
        step_name_pattern: str
            Pattern for multiple inlists (e.g., 'inlist_project_step{}')
        """
        if not inlists:
            return

        os.makedirs(output_path, exist_ok=True)

        if len(inlists) == 1:
            inlists[0].to_file(f'{output_path}/{single_name}')
        else:
            for i, inlist in enumerate(inlists):
                inlist.to_file(f'{output_path}/{step_name_pattern.format(i)}')

    def write_inlists(self, output_dir):
        """Write all managed inlists to the output directory.

        Parameters
        ----------
        output_dir : str
            Path to the output directory.
        """
        self._write_inlist_group(self.binary_inlists, f'{output_dir}/binary',
                                  'inlist_project', 'inlist_project_step{}')
        self._write_inlist_group(self.binary_star1_inlists, f'{output_dir}/binary',
                                  'inlist1', 'inlist1_step{}')
        self._write_inlist_group(self.binary_star2_inlists, f'{output_dir}/binary',
                                  'inlist2', 'inlist2_step{}')
        self._write_inlist_group(self.star1_inlists, f'{output_dir}/star1',
                                  'inlist_step0', 'inlist_step{}')
        self._write_inlist_group(self.star2_inlists, f'{output_dir}/star2',
                                  'inlist_step0', 'inlist_step{}')
