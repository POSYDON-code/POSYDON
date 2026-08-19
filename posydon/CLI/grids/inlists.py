"""Inlist abstraction and inlist builder functions for the grid setup."""

import copy
import os
from dataclasses import dataclass

from posydon.utils.posydonwarning import Pwarn


class InlistSection:
    """Represents a single namelist section (e.g., &star_job, &controls)."""

    def __init__(self, name: str, parameters: dict = None):
        self.name = name
        self.parameters = parameters if parameters is not None else {}

    def __repr__(self):
        """Return a string representation of the section."""
        param_count = len(self.parameters)
        if param_count == 0:
            return f"InlistSection(name='{self.name}', parameters=0)"

        preview_keys = list(self.parameters.keys())[:3]
        preview = ", ".join(preview_keys)
        if param_count > 3:
            preview += ", ..."

        return f"InlistSection(name='{self.name}', parameters={param_count}: [{preview}])"

    def merge(self, other: "InlistSection") -> "InlistSection":
        """Merge another section into this one (later values override)."""
        if self.name != other.name:
            raise ValueError("Cannot merge sections with different names")
        self.parameters.update(other.parameters)
        return self

    def to_fortran(self) -> str:
        """Convert the section to a Fortran namelist format string."""
        lines = [f"&{self.name}"]
        for key, value in self.parameters.items():
            lines.append(f"\t{key} = {self._format_value(value)}")
        lines.append(f"/ ! end of {self.name} namelist\n")
        return "\n".join(lines)

    @staticmethod
    def _format_value(value) -> str:
        """Format a value according to the legacy run-grid Fortran rule."""
        if isinstance(value, bool):
            return '.true.' if value else '.false.'
        elif isinstance(value, float):
            if value.is_integer():
                return str(int(value))
            return '{0:.10e}'.format(value).replace('e', 'd')
        return value

    @staticmethod
    def _parse_value(value_str: str) -> str:
        """Parse a value string by stripping surrounding whitespace."""
        return value_str.strip()

    @classmethod
    def from_string(cls, name: str, content: str) -> "InlistSection":
        """Parse a namelist section from string content."""
        parameters = {}

        for line in content.split('\n'):
            if '!' in line:
                line = line[:line.index('!')]
            line = line.strip()

            if not line or line.startswith('&') or line == '/':
                continue

            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                parameters[key] = cls._parse_value(value)

        return cls(name, parameters)

    @classmethod
    def from_file(cls, name: str, path: str) -> "InlistSection":
        """Parse a namelist section from a file."""
        with open(path, 'r') as f:
            content = f.read()
        return cls.from_string(name, content)


class Inlist:
    """Represents a complete inlist file with multiple sections."""

    def __init__(self, name: str, sections: dict = None):
        object.__setattr__(self, 'name', name)
        object.__setattr__(self, 'sections', sections if sections is not None else {})

    def __repr__(self):
        """Return a string representation of the inlist."""
        section_count = len(self.sections)
        if section_count == 0:
            return f"Inlist(name='{self.name}', sections=0)"

        section_names = list(self.sections.keys())
        sections_str = ", ".join(section_names)

        return f"Inlist(name='{self.name}', sections={section_count}: [{sections_str}])"

    def __getattr__(self, name):
        """Allow attribute-style access to sections (e.g., inlist.controls)."""
        sections = object.__getattribute__(self, 'sections')
        if name in sections:
            return sections[name]

        raise AttributeError(f"Inlist has no section '{name}'")

    def __setattr__(self, name, value):
        """Allow setting sections as attributes."""
        if name in ['name', 'sections']:
            object.__setattr__(self, name, value)
        elif isinstance(value, InlistSection):
            self.sections[name] = value
        else:
            object.__setattr__(self, name, value)

    def add_section(self, section: InlistSection) -> None:
        """Add or merge a section into the inlist."""
        if section.name in self.sections:
            self.sections[section.name].merge(section)
        else:
            self.sections[section.name] = section

    def merge(self, other: "Inlist") -> "Inlist":
        """Merge another inlist into a new inlist (later values override)."""
        merged = Inlist(self.name, copy.deepcopy(self.sections))
        for section in other.sections.values():
            merged.add_section(section)
        return merged

    def to_file(self, filepath: str) -> None:
        """Write the inlist to a file with a fixed section ordering."""
        lines = []

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
    def from_file(cls, filepath: str, name: str = None, section: str = None) -> "Inlist":
        """Read and parse an inlist file.

        Args:
            filepath: Path to the inlist file.
            name: Optional name for the inlist (defaults to the filename).
            section: Optional section name for files without section markers.

        Returns:
            An Inlist object with the parsed sections.
        """
        if name is None:
            name = os.path.basename(filepath)

        with open(filepath, 'r') as f:
            content = f.read()

        return cls.from_string(content, name, section=section)

    @classmethod
    def from_string(
        cls,
        content: str,
        name: str = "inlist",
        section: str = None,
    ) -> "Inlist":
        """Parse an inlist from string content following clean_inlist_file.

        Args:
            content: String content of the inlist file.
            name: Name for the inlist.
            section: Optional section name if the file has no section markers.

        Returns:
            An Inlist object with the parsed sections.
        """
        inlist = cls(name)

        param_value_list = []
        for line in content.split('\n'):
            param_and_value = line.strip('\n').strip().split('!')[0].strip()
            if param_and_value and ('=' in param_and_value or '&' in param_and_value):
                param_value_list.append(param_and_value)

        sections_found = {k: {} for k in param_value_list if '&' in k}

        if not sections_found:
            if section:
                parameters = {}
                for item in param_value_list:
                    key, value = item.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    if value not in ["''", "'.'"]:
                        parameters[key] = InlistSection._parse_value(value)
                inlist.add_section(InlistSection(section, parameters))
        else:
            current_section = None
            for item in param_value_list:
                if '&' in item:
                    current_section = item.replace('&', '').strip()
                    sections_found[item] = {}
                elif current_section:
                    key, value = item.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    if value not in ["''", "'.'"]:
                        sections_found['&' + current_section][key] = InlistSection._parse_value(value)

            for section_key, params in sections_found.items():
                section_name = section_key.replace('&', '').strip()
                inlist.add_section(InlistSection(section_name, params))

        return inlist


class MESAInlists:
    """Handles MESA inlists for single and binary star evolution."""

    def __init__(self, path: str):
        self.path = path

        controls_inlist = Inlist.from_file(f'{self.path}/star/controls.defaults', section='controls')
        star_job_inlist = Inlist.from_file(f'{self.path}/star/star_job.defaults', section='star_job')
        self.base_star_inlist = controls_inlist.merge(star_job_inlist)

        binary_job_inlist = Inlist.from_file(f'{self.path}/binary/binary_job.defaults', section='binary_job')
        binary_controls_inlist = Inlist.from_file(f'{self.path}/binary/binary_controls.defaults', section='binary_controls')
        self.base_binary_inlist = binary_job_inlist.merge(binary_controls_inlist)

        sections_to_clean = [
            self.base_star_inlist.controls,
            self.base_star_inlist.star_job,
            self.base_binary_inlist.binary_controls,
            self.base_binary_inlist.binary_job,
        ]
        for section in sections_to_clean:
            section.parameters = self._clean_parameters(section.parameters)

    def _clean_parameters(self, params_dict: dict) -> dict:
        """Clean parameters by removing read_extra/inlist refs and renaming num_x_ctrls.

        Args:
            params_dict: Dictionary of parameters to clean.

        Returns:
            The cleaned parameters dictionary.
        """
        cleaned = {k: v for k, v in params_dict.items()
                   if not any(substring in k for substring in ['read_extra', 'inlist'])}

        keys_to_replace = {k: k.replace('num_x_ctrls', '1')
                           for k in cleaned.keys() if 'num_x_ctrls' in k}
        for old_key, new_key in keys_to_replace.items():
            cleaned[new_key] = cleaned.pop(old_key)

        return cleaned

    def __repr__(self):
        """Return a string representation of the MESA inlists."""
        return (f"MESAInlists(path='{self.path}', "
                f"base_star_inlist={self.base_star_inlist}, "
                f"base_binary_inlist={self.base_binary_inlist})")


class InlistManager:
    """Manages inlists for different evolution steps and phases in MESA."""

    def __init__(self):
        self.binary_inlists = []
        self.binary_star1_inlists = []
        self.binary_star2_inlists = []
        self.star1_inlists = []
        self.star2_inlists = []

    def keys(self) -> list:
        """Return the names of the managed inlist groups."""
        return ['binary_inlists', 'binary_star1_inlists', 'binary_star2_inlists', 'star1_inlists', 'star2_inlists']

    def __getitem__(self, key):
        """Return the inlist group for the given key."""
        if key in self.keys():
            return getattr(self, key)
        else:
            raise KeyError(f"InlistManager has no key '{key}'")

    def append_binary_inlist(self, inlist: Inlist) -> None:
        """Append an inlist to the binary_inlists group."""
        self.binary_inlists.append(inlist)

    def append_binary_star1_inlist(self, inlist: Inlist) -> None:
        """Append an inlist to the binary_star1_inlists group."""
        self.binary_star1_inlists.append(inlist)

    def append_binary_star2_inlist(self, inlist: Inlist) -> None:
        """Append an inlist to the binary_star2_inlists group."""
        self.binary_star2_inlists.append(inlist)

    def append_star1_inlist(self, inlist: Inlist) -> None:
        """Append an inlist to the star1_inlists group."""
        self.star1_inlists.append(inlist)

    def append_star2_inlist(self, inlist: Inlist) -> None:
        """Append an inlist to the star2_inlists group."""
        self.star2_inlists.append(inlist)

    def __repr__(self):
        """Return a string representation of the inlist manager."""
        return (f"InlistManager(\nbinary_inlists={len(self.binary_inlists)}, \n"
                f"binary_star1_inlists={len(self.binary_star1_inlists)}, \n"
                f"binary_star2_inlists={len(self.binary_star2_inlists)}, \n"
                f"star1_inlists={len(self.star1_inlists)}, \n"
                f"star2_inlists={len(self.star2_inlists)})")

    def _write_inlist_group(
        self,
        inlists: list,
        output_path: str,
        single_name: str,
        step_name_pattern: str,
    ) -> None:
        """Write a group of inlists with consistent naming logic.

        Args:
            inlists: List of Inlist objects to write.
            output_path: Directory path to write to.
            single_name: Filename if only one inlist (e.g., 'inlist_project').
            step_name_pattern: Pattern for multiple inlists (e.g., 'inlist_step{}').
        """
        if not inlists:
            return

        os.makedirs(output_path, exist_ok=True)

        if len(inlists) == 1:
            inlists[0].to_file(f'{output_path}/{single_name}')
        else:
            for i, inlist in enumerate(inlists):
                inlist.to_file(f'{output_path}/{step_name_pattern.format(i)}')

    def write_inlists(self, output_dir: str) -> None:
        """Write all managed inlists to the output directory.

        Args:
            output_dir: Path to the output directory.
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


@dataclass
class MergedInlists:
    """Dataclass holding the merged inlist dictionaries and grid parameter lists."""

    final_binary_controls: dict
    final_binary_job: dict
    final_star1_binary_controls: dict
    final_star1_binary_job: dict
    final_star2_binary_controls: dict
    final_star2_binary_job: dict
    grid_params_binary_controls: list
    grid_params_binary_job: list
    grid_params_star1_binary_controls: list
    grid_params_star1_binary_job: list
    grid_params_star2_binary_controls: list
    grid_params_star2_binary_job: list


def merge_inlist_layers(
    mesa_inlists: dict,
    grid_parameters: list,
    working_directory: str,
) -> MergedInlists:
    """Merge the various inlist layers for a binary grid.

    Args:
        mesa_inlists: All of the values from the mesa_inlists section of the
            inifile.
        grid_parameters: A list of the parameters from the csv file.
        working_directory: The working directory where the grid is being set up.

    Returns:
        MergedInlists with the final inlist dictionaries and grid parameter
        lists.
    """
    if mesa_inlists['single_star_grid']:
        return MergedInlists({}, {}, {}, {}, {}, {}, [], [], [], [], [], [])

    final_binary_controls = {}
    final_binary_job = {}
    final_star1_binary_job = {}
    final_star2_binary_job = {}
    final_star1_binary_controls = {}
    final_star2_binary_controls = {}

    for k, v in mesa_inlists.items():
        if v is not None:
            if 'binary_controls' in k:
                controls_dict = Inlist.from_file(v, section='binary_controls')
                for k1, v1 in controls_dict.binary_controls.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    final_binary_controls[k1] = v1

            elif 'binary_job' in k:
                job_dict = Inlist.from_file(v, section='binary_job')
                for k1, v1 in job_dict.binary_job.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    final_binary_job[k1] = v1

            elif 'star1_job' in k:
                star_job1_dict = Inlist.from_file(v, section='star_job')
                for k1, v1 in star_job1_dict.star_job.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    final_star1_binary_job[k1] = v1

            elif 'star2_job' in k:
                star_job2_dict = Inlist.from_file(v, section='star_job')
                for k1, v1 in star_job2_dict.star_job.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    final_star2_binary_job[k1] = v1

            elif 'star1_controls' in k:
                star_control1_dict = Inlist.from_file(v, section='controls')
                for k1, v1 in star_control1_dict.controls.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    if 'num_x_ctrls' in k1:
                        final_star1_binary_controls[k1.replace('num_x_ctrls', '1')] = v1
                    else:
                        final_star1_binary_controls[k1] = v1

            elif 'star2_controls' in k:
                star_control2_dict = Inlist.from_file(v, section='controls')
                for k1, v1 in star_control2_dict.controls.parameters.items():
                    if ('read_extra' in k1) or ('inlist' in k1):
                        continue
                    if 'num_x_ctrls' in k1:
                        final_star2_binary_controls[k1.replace('num_x_ctrls', '1')] = v1
                    else:
                        final_star2_binary_controls[k1] = v1

    grid_params_binary_controls = [param for param in grid_parameters if param in final_binary_controls.keys()]
    print("Grid parameters that effect binary_controls: {0}".format(','.join(grid_params_binary_controls)))
    grid_params_binary_job = [param for param in grid_parameters if param in final_binary_job.keys()]
    print("Grid parameters that effect binary_job: {0}".format(','.join(grid_params_binary_job)))

    grid_params_star1_binary_controls = [param for param in grid_parameters if param in final_star1_binary_controls.keys()]
    print("Grid parameters that effect star1_binary_controls: {0}".format(','.join(grid_params_star1_binary_controls)))
    grid_params_star1_binary_job = [param for param in grid_parameters if param in final_star1_binary_job.keys()]
    print("Grid parameters that effect star1_binary_job: {0}".format(','.join(grid_params_star1_binary_job)))

    grid_params_star2_binary_controls = [param for param in grid_parameters if param in final_star2_binary_controls.keys()]
    print("Grid parameters that effect star2_binary_controls: {0}".format(','.join(grid_params_star2_binary_controls)))
    grid_params_star2_binary_job = [param for param in grid_parameters if param in final_star2_binary_job.keys()]
    print("Grid parameters that effect star2_binary_job: {0}".format(','.join(grid_params_star2_binary_job)))

    if grid_params_star1_binary_controls:
        final_star1_binary_controls['read_extra_controls_inlist1'] = '.true.'
        final_star1_binary_controls['extra_controls_inlist1_name'] = "'inlist_grid_star1_binary_controls'"

    if grid_params_star2_binary_controls:
        final_star2_binary_controls['read_extra_controls_inlist1'] = '.true.'
        final_star2_binary_controls['extra_controls_inlist1_name'] = "'inlist_grid_star2_binary_controls'"

    if grid_params_star1_binary_job:
        final_star1_binary_job['read_extra_star_job_inlist1'] = '.true.'
        final_star1_binary_job['extra_star_job_inlist1_name'] = "'inlist_grid_star1_binary_job'"

    if grid_params_star2_binary_job:
        final_star2_binary_job['read_extra_star_job_inlist1'] = '.true.'
        final_star2_binary_job['extra_star_job_inlist1_name'] = "'inlist_grid_star2_binary_job'"

    final_binary_job['inlist_names(1)'] = "'{0}'".format(os.path.join(working_directory, 'binary', 'inlist1'))
    final_binary_job['inlist_names(2)'] = "'{0}'".format(os.path.join(working_directory, 'binary', 'inlist2'))

    return MergedInlists(
        final_binary_controls,
        final_binary_job,
        final_star1_binary_controls,
        final_star1_binary_job,
        final_star2_binary_controls,
        final_star2_binary_job,
        grid_params_binary_controls,
        grid_params_binary_job,
        grid_params_star1_binary_controls,
        grid_params_star1_binary_job,
        grid_params_star2_binary_controls,
        grid_params_star2_binary_job,
    )


def build_star_formation(
    star_key: str,
    mesa_inlists: dict,
    working_directory: str,
    counterpart_binary_job: dict,
) -> str | None:
    """Construct the per-step formation inlists for a star.

    Args:
        star_key: Either 'star1' or 'star2'.
        mesa_inlists: All of the values from the mesa_inlists section of the
            inifile.
        working_directory: The working directory where the grid is being set up.
        counterpart_binary_job: The final binary job dictionary for this star,
            which is wired to load the model produced by the last step.

    Returns:
        The space-joined list of written step inlist paths, or None when no
        formation inlists were supplied.
    """
    star_formation = {}

    if star_key == 'star1':
        if (mesa_inlists['zams_filename_1'] is not None) and (not mesa_inlists['single_star_grid']):
            star_formation_dictionary = {}
        elif mesa_inlists['single_star_grid']:
            star_formation_dictionary = dict(filter(lambda elem: (('star1_job' in elem[0]) or ('star1_controls' in elem[0])) and elem[1] is not None, mesa_inlists.items()))
        else:
            star_formation_dictionary = dict(filter(lambda elem: 'star1_formation' in elem[0] and elem[1] is not None, mesa_inlists.items()))
    else:
        if mesa_inlists['zams_filename_2'] is not None:
            star_formation_dictionary = {}
        else:
            star_formation_dictionary = dict(filter(lambda elem: 'star2_formation' in elem[0] and elem[1] is not None, mesa_inlists.items()))

    if star_formation_dictionary:
        inlist_star_formation = ''
        number_of_steps = 1
        for k, v in star_formation_dictionary.items():
            if type(v) == list:
                number_of_steps = max(number_of_steps, len(v))

        for step in range(number_of_steps):
            star_formation['step{0}'.format(step)] = {}
            star_formation['step{0}'.format(step)]['inlist_file'] = os.path.join(working_directory, star_key, 'inlist_step{0}'.format(step))
            for k, v in star_formation_dictionary.items():
                star_formation['step{0}'.format(step)][k] = v[step] if type(v) == list else v

        for step, step_inlists in enumerate(star_formation.values()):
            final_controls = {}
            final_job = {}
            for k, v in step_inlists.items():
                if star_key == 'star1':
                    controls_match = ('star1_formation_controls' in k) or ('star1_controls' in k)
                    job_match = ('star1_formation_job' in k) or ('star1_job' in k)
                else:
                    controls_match = 'star2_formation_controls' in k
                    job_match = 'star2_formation_job' in k
                if controls_match:
                    controls_dict = Inlist.from_file(v, section='controls')
                    for k1, v1 in controls_dict.controls.parameters.items():
                        if ('read_extra' in k1) or ('inlist' in k1):
                            continue
                        if 'num_x_ctrls' in k1:
                            final_controls[k1.replace('num_x_ctrls', '1')] = v1
                        else:
                            final_controls[k1] = v1

                if job_match:
                    job_dict = Inlist.from_file(v, section='star_job')
                    for k1, v1 in job_dict.star_job.parameters.items():
                        if ('read_extra' in k1) or ('inlist' in k1):
                            continue
                        final_job[k1] = v1

            counterpart_binary_job['create_pre_main_sequence_model'] = ".false."
            counterpart_binary_job['load_saved_model'] = ".true."
            counterpart_binary_job['saved_model_name'] = "'initial_{0}_step{1}.mod'".format(star_key, step)
            if step == 0:
                final_job['save_model_when_terminate'] = '.true.'
                final_job['save_model_filename'] = "'initial_{0}_step{1}.mod'".format(star_key, step)
            else:
                final_job['create_pre_main_sequence_model'] = ".false."
                final_job['load_saved_model'] = ".true."
                final_job['saved_model_name'] = "'initial_{0}_step{1}.mod'".format(star_key, step - 1)
                final_job['save_model_when_terminate'] = '.true.'
                final_job['save_model_filename'] = "'initial_{0}_step{1}.mod'".format(star_key, step)

            if star_key == 'star1':
                if (mesa_inlists['zams_filename_1'] is not None) and (mesa_inlists['single_star_grid']):
                    final_controls['zams_filename'] = "'{0}'".format(mesa_inlists['zams_filename_1'])
                elif (mesa_inlists['zams_filename_1'] is None) and (mesa_inlists['single_star_grid']) and ('zams_filename_1' in final_controls.keys()):
                    final_controls.pop("zams_filename", None)

            if os.path.exists(step_inlists['inlist_file']):
                Pwarn('Replace ' + step_inlists['inlist_file'], "OverwriteWarning")
            with open(step_inlists['inlist_file'], 'wb') as f:
                f.write(b'&controls\n\n')
                for k, v in final_controls.items():
                    f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

                f.write(b'\n\n')

                f.write("""
        / ! end of {0}_controls namelist

        """.format(star_key).encode('utf-8'))

                f.write(b'&star_job\n\n')
                for k, v in final_job.items():
                    f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

                f.write("""
        / ! end of {0}_job namelist

        """.format(star_key).encode('utf-8'))

            inlist_star_formation += ' {0}'.format(step_inlists['inlist_file'])

        return inlist_star_formation
    else:
        return None


def apply_output_controls(
    final_binary_controls: dict,
    final_star1_binary_controls: dict,
    final_star2_binary_controls: dict,
    final_star1_binary_job: dict,
    final_star2_binary_job: dict,
    mesa_inlists: dict,
) -> None:
    """Apply the MESA binary output control flags to the final inlist dicts.

    Args:
        final_binary_controls: The merged binary controls dictionary.
        final_star1_binary_controls: The merged star1 controls dictionary.
        final_star2_binary_controls: The merged star2 controls dictionary.
        final_star1_binary_job: The merged star1 job dictionary.
        final_star2_binary_job: The merged star2 job dictionary.
        mesa_inlists: All of the values from the mesa_inlists section of the
            inifile.
    """
    if mesa_inlists['final_profile_star1']:
        final_star1_binary_job['write_profile_when_terminate'] = ".true."
        final_star1_binary_job['filename_for_profile_when_terminate'] = "'final_profile_star1.data'"
    else:
        final_star1_binary_job['write_profile_when_terminate'] = ".false."

    if mesa_inlists['final_profile_star2']:
        final_star2_binary_job['write_profile_when_terminate'] = ".true."
        final_star2_binary_job['filename_for_profile_when_terminate'] = "'final_profile_star2.data'"
    else:
        final_star2_binary_job['write_profile_when_terminate'] = ".false."

    if mesa_inlists['final_model_star1']:
        final_star1_binary_job['save_model_when_terminate'] = ".true."
        final_star1_binary_job['save_model_filename'] = "'final_star1.mod'"
    else:
        final_star1_binary_job['save_model_when_terminate'] = ".false."

    if mesa_inlists['final_model_star2']:
        final_star2_binary_job['save_model_when_terminate'] = ".true."
        final_star2_binary_job['save_model_filename'] = "'final_star2.mod'"
    else:
        final_star2_binary_job['save_model_when_terminate'] = ".false."

    if mesa_inlists['history_star1']:
        final_star1_binary_controls['do_history_file'] = ".true."
    else:
        final_star1_binary_controls['do_history_file'] = ".false."

    if mesa_inlists['history_star2']:
        final_star2_binary_controls['do_history_file'] = ".true."
    else:
        final_star2_binary_controls['do_history_file'] = ".false."

    final_binary_controls['history_interval'] = mesa_inlists['history_interval']
    final_star1_binary_controls['history_interval'] = mesa_inlists['history_interval']
    final_star2_binary_controls['history_interval'] = mesa_inlists['history_interval']

    if not mesa_inlists['binary_history']:
        final_binary_controls['history_interval'] = "-1"

    if mesa_inlists['zams_filename_1'] is not None:
        final_star1_binary_controls['zams_filename'] = "'{0}'".format(mesa_inlists['zams_filename_1'])
    if mesa_inlists['zams_filename_2'] is not None:
        final_star2_binary_controls['zams_filename'] = "'{0}'".format(mesa_inlists['zams_filename_2'])


def write_binary_inlists(
    working_directory: str,
    final_binary_controls: dict,
    final_binary_job: dict,
    final_star1_binary_controls: dict,
    final_star1_binary_job: dict,
    final_star2_binary_controls: dict,
    final_star2_binary_job: dict,
) -> tuple[str, str, str]:
    """Write the binary inlist_project, inlist1 and inlist2 files.

    Args:
        working_directory: The working directory where the binary directory
            lives.
        final_binary_controls: The merged binary controls dictionary.
        final_binary_job: The merged binary job dictionary.
        final_star1_binary_controls: The merged star1 controls dictionary.
        final_star1_binary_job: The merged star1 job dictionary.
        final_star2_binary_controls: The merged star2 controls dictionary.
        final_star2_binary_job: The merged star2 job dictionary.

    Returns:
        The paths to inlist_project, inlist1 and inlist2.
    """
    inlist_binary_project = os.path.join(working_directory, 'binary', 'inlist_project')
    inlist_star1_binary = os.path.join(working_directory, 'binary', 'inlist1')
    inlist_star2_binary = os.path.join(working_directory, 'binary', 'inlist2')

    if os.path.exists(inlist_binary_project):
        Pwarn('Replace ' + inlist_binary_project, "OverwriteWarning")
    with open(inlist_binary_project, 'wb') as f:
        f.write(b'&binary_controls\n\n')
        for k, v in final_binary_controls.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b'\n/ ! end of binary_controls namelist')

        f.write(b'\n\n')

        f.write(b'&binary_job\n\n')
        for k, v in final_binary_job.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b"""
/ ! end of binary_job namelist
                """)

    if os.path.exists(inlist_star1_binary):
        Pwarn('Replace ' + inlist_star1_binary, "OverwriteWarning")
    with open(inlist_star1_binary, 'wb') as f:
        f.write(b'&controls\n\n')
        for k, v in final_star1_binary_controls.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b'\n\n')

        f.write(b"""
/ ! end of star1_controls namelist

""")

        f.write(b'&star_job\n\n')
        for k, v in final_star1_binary_job.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b"""
/ ! end of star1_job namelist

""")

    if os.path.exists(inlist_star2_binary):
        Pwarn('Replace ' + inlist_star2_binary, "OverwriteWarning")
    with open(inlist_star2_binary, 'wb') as f:
        f.write(b'&controls\n\n')
        for k, v in final_star2_binary_controls.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b'\n\n')

        f.write(b"""
/ ! end of star2_controls namelist

""")

        f.write(b'&star_job\n\n')
        for k, v in final_star2_binary_job.items():
            f.write('\t{0} = {1}\n'.format(k, v).encode('utf-8'))

        f.write(b"""
/ ! end of star2_job namelist

""")

    return inlist_binary_project, inlist_star1_binary, inlist_star2_binary


def construct_static_inlist(
    mesa_inlists: dict,
    grid_parameters: list,
    working_directory: str = os.getcwd(),
) -> tuple[str, str, str, str, str]:
    """Based on all the inlists that were passed construct the MESA project dir.

    Args:
        mesa_inlists: All of the values from the mesa_inlists section of the
            inifile.
        grid_parameters: A list of the parameters from the csv file.
        working_directory: The working directory where the grid is being set up.

    Returns:
        (inlist_star1_formation, inlist_star2_formation, inlist_binary_project,
         inlist_star1_binary, inlist_star2_binary).
    """
    if 'zams_filename' in mesa_inlists.keys():
        mesa_inlists['zams_filename_1'] = mesa_inlists['zams_filename']
        mesa_inlists['zams_filename_2'] = mesa_inlists['zams_filename']

    if 'zams_filename_1' not in mesa_inlists.keys():
        mesa_inlists['zams_filename_1'] = None

    if 'zams_filename_2' not in mesa_inlists.keys():
        mesa_inlists['zams_filename_2'] = None

    if 'single_star_grid' not in mesa_inlists.keys():
        mesa_inlists['single_star_grid'] = False

    merged = merge_inlist_layers(mesa_inlists, grid_parameters, working_directory)

    inlist_star1_formation = build_star_formation(
        'star1', mesa_inlists, working_directory, merged.final_star1_binary_job
    )
    inlist_star2_formation = build_star_formation(
        'star2', mesa_inlists, working_directory, merged.final_star2_binary_job
    )

    apply_output_controls(
        merged.final_binary_controls,
        merged.final_star1_binary_controls,
        merged.final_star2_binary_controls,
        merged.final_star1_binary_job,
        merged.final_star2_binary_job,
        mesa_inlists,
    )

    inlist_binary_project, inlist_star1_binary, inlist_star2_binary = write_binary_inlists(
        working_directory,
        merged.final_binary_controls,
        merged.final_binary_job,
        merged.final_star1_binary_controls,
        merged.final_star1_binary_job,
        merged.final_star2_binary_controls,
        merged.final_star2_binary_job,
    )

    return (
        inlist_star1_formation,
        inlist_star2_formation,
        inlist_binary_project,
        inlist_star1_binary,
        inlist_star2_binary,
    )
