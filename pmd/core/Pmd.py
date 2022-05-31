import inspect
import os
from typing import Any, Dict, List, Optional, Union

import yaml

import pmd.core
# These have to be written explicitly for typing
from pmd.core.ForceField import ForceField
from pmd.core.Job import Job
from pmd.core.Lammps import Lammps
from pmd.core.Procedure import Procedure
from pmd.core.System import System
from pmd.util import Pmdlogging, build_dir

PRIMITIVE_TYPES = (str, int, float, bool)
SUPPORTED_YAML_EXTS = (".yml", ".yaml")
OBJECT_PRFIX = 'pmd.'


# Custom yaml config file dictionary constructor
def to_yaml_dict(
        cls: Union[System, ForceField, Lammps, Procedure, Job]) -> Dict:
    return {
        # strip off the front underscore and only add to dict
        # if value is not None
        k.lstrip('_'): custom_class_yaml_dumper(v)
        for k, v in cls.__dict__.items() if v is not None
    }


# Custom method to dump values of non-primitive type to the yaml config file
def custom_class_yaml_dumper(v: Any) -> Any:
    return_value = v

    # If value is a list, recursively go through each item in the list
    # Specifically, this is for the Lammps procedure list
    if isinstance(v, list):
        return_value = [custom_class_yaml_dumper(i) for i in v]

    # If value is a non-primitive type, expand it to a dict
    elif not isinstance(v, PRIMITIVE_TYPES):
        return_value = {f"{OBJECT_PRFIX}{v}": to_yaml_dict(v)}

    return return_value


def instantiate_from_cls_name(class_name: str, prop_dict: dict):
    # first obtain a list of all classes in this module
    class_list = inspect.getmembers(pmd.core, inspect.isclass)
    class_dict = {k: v for k, v in class_list}

    # find the matching class
    filtered_class_name = class_name.lstrip(OBJECT_PRFIX).split('-')[0]
    the_class = class_dict.get(filtered_class_name, None)
    if the_class is None:
        raise NameError(
            f'{class_name} type is not found in {pmd.core.__name__} module')

    # get the constructor parameter list of the class
    sig = inspect.signature(the_class.__init__)
    param_keys = list(sig.parameters.keys())
    if param_keys[0] == 'self':
        param_keys = param_keys[1:]

    # remove props not in the parameter list of the class
    filtered_prop_dict = {
        k: custom_class_yaml_loader(v)
        for k, v in prop_dict.items() if k in param_keys
    }

    Pmdlogging.info(
        f'{class_name} object successfully loaded from the YAML file.')
    return the_class(**filtered_prop_dict)


# Custom method to load values from the yaml config file
def custom_class_yaml_loader(v: Any) -> Any:
    return_value = v

    # If value is a list, recursively go through each item in the list
    # Specifically, this is for the Lammps procedure list
    if isinstance(v, list):
        return_value = [custom_class_yaml_loader(i) for i in v]

    # If value is a dict, instantiate it to an object
    elif isinstance(v, dict):
        class_name, props_dict = next(iter(v.items()))
        return_value = instantiate_from_cls_name(class_name, props_dict)

    # If value is starts with pmd., instantiate it to an object with
    # default params
    elif isinstance(v, str) and v.startswith(OBJECT_PRFIX):
        return_value = instantiate_from_cls_name(v, {})

    return return_value


class Pmd:
    '''Template object to perform tasks for Systems, Lammps, and Jobs
       altogether (e.g. create data files, lammps input files, job scheduler
       files, or config files)

        Attributes:
            system (System): a System object

            lammps (Lammps or list[Lammps]): one or a list of Lammps objects

            job (Job or list[Job]): one or a list of Job objects
    '''

    def __init__(
        self,
        system: Optional[System] = None,
        lammps: Optional[Union[Lammps, List[Lammps]]] = None,
        job: Optional[Union[Job, List[Job]]] = None,
    ):
        if lammps and not isinstance(lammps, list):
            lammps = [lammps]
        if job and not isinstance(job, list):
            job = [job]

        self._system = system
        self._lammps = lammps
        self._job = job

    @build_dir
    def create(self,
               output_dir: str = '.',
               save_config: bool = False,
               config_fname: str = 'config.yaml') -> None:
        '''Method to create files from all the pmd objects. This method can
        can also automatically generate a config file if `save_config` input
        argument is set to True.

        Parameters:
            output_dir (str): Directory for all the generated files; default:
                              `"."`

            save_config (bool): Whether to save a config file; default: `False`

            config_fname (str): Name of the config file; default:
                                `"config.yaml"`

        Returns:
            None
        '''
        if self._system:
            self._system.write_data(output_dir)
        if self._lammps:
            for lmp in self._lammps:
                lmp.write_lammps(output_dir)
        if self._job:
            for job in self._job:
                job.write_job(output_dir)

        if save_config:
            self.save_config(output_dir, config_fname)

    @build_dir
    def save_config(self, output_dir: str, config_fname: str = 'config.yaml'):
        '''Method to create a config file with all the details of the System,
        Lammps, or Job settings. This method only creates the config file.

        Parameters:
            output_dir (str): Directory for all the generated files; default:
                              `"."`

            config_fname (str): Name of the config file; default:
                                `"config.yaml"`

        Returns:
            None
        '''

        config_dict = {}

        if self._system:
            config_dict[f'{OBJECT_PRFIX}{self._system}'] = to_yaml_dict(
                self._system)

        if self._lammps and len(self._lammps) > 1:
            for i, lmp in enumerate(self._lammps):
                config_dict[f'{OBJECT_PRFIX}{lmp}-{i}'] = to_yaml_dict(lmp)
        elif self._lammps and len(self._lammps) == 1:
            config_dict[f'{OBJECT_PRFIX}{self._lammps[0]}'] = to_yaml_dict(
                self._lammps[0])

        if self._job and len(self._job) > 1:
            for i, job in enumerate(self._job):
                config_dict[f'{OBJECT_PRFIX}{job}-{i}'] = to_yaml_dict(job)
        elif self._job and len(self._job) == 1:
            config_dict[f'{OBJECT_PRFIX}{self._job[0]}'] = to_yaml_dict(
                self._job[0])

        with open(os.path.join(output_dir, config_fname), 'w') as yaml_file:
            yaml.safe_dump(config_dict, yaml_file, sort_keys=False)

        Pmdlogging.info(f'Config file - {config_fname} successfully '
                        f'saved to {output_dir}')

    @staticmethod
    def load_config(config_file: str, output_dir: str = '.'):
        '''Method to load a config file and create all the objects listed in
        the config file

        Parameters:
            config_file (str): Config file to load

            output_dir (str): Directory for all the generated files; default:
                              `"."`

        Returns:
            None
        '''

        if os.path.splitext(config_file)[1] not in SUPPORTED_YAML_EXTS:
            raise ValueError(
                f'The file you are loading does not seem to be a yaml file'
                f'(file must end with {" ,".join(SUPPORTED_YAML_EXTS)})')

        with open(config_file) as yaml_file:
            yaml_dict = yaml.safe_load(yaml_file)
            for k, v in yaml_dict.items():
                obj = instantiate_from_cls_name(k, v)
                if isinstance(obj, System):
                    obj.write_data(output_dir)
                elif isinstance(obj, Lammps):
                    obj.write_lammps(output_dir)
                elif isinstance(obj, Job):
                    obj.write_job(output_dir)
