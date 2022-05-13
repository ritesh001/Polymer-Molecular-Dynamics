import os
import yaml
from typing import Any, Dict, List, Union

from pmd.core.ForceField import ForceField
from pmd.core.Job import Job
from pmd.core.Lammps import Lammps
from pmd.core.Procedure import Procedure
from pmd.core.System import System
from pmd.util import Util

PRIMITIVE_TYPES = (str, int, float, bool)
EXCLUDED_VALUE_TYPES = (System, Lammps)
SUPPORTED_YAML_EXTS = (".yml", ".yaml")


# Custom yaml config file dictionary constructor
def to_yaml_dict(
        cls: Union[System, ForceField, Lammps, Procedure, Job]) -> Dict:
    return {
        # strip off the front underscore
        k.lstrip('_'): custom_class_yaml_dumper(v)
        for k, v in cls.__dict__.items()
        # only add to dict if value is not None nor System/Lammps instances
        if v is not None and not issubclass(type(v), EXCLUDED_VALUE_TYPES)
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
        return_value = {str(v): to_yaml_dict(v)}
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
        system: System = None,
        lammps: Union[Lammps, List[Lammps]] = None,
        job: Union[Job, List[Job]] = None,
    ):
        if lammps and not isinstance(lammps, list):
            lammps = [lammps]
        if job and not isinstance(job, list):
            job = [job]

        self._system = system
        self._lammps = lammps
        self._job = job
        self._config_items: List[Union[System, Lammps, Job]] = []

    @Util.build_dir
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
            self._config_items.append(self._system)
        if self._lammps:
            for lmp in self._lammps:
                lmp.write_lammps(output_dir)
                self._config_items.append(lmp)
        if self._job:
            for job in self._job:
                job.write_job(output_dir)
                self._config_items.append(job)

        if save_config:
            self.save_config(output_dir, config_fname)

    @Util.build_dir
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
        with open(os.path.join(output_dir, config_fname), 'w') as yaml_file:
            config = {
                str(config_item): to_yaml_dict(config_item)
                for config_item in self._config_items
            }
            yaml.dump(config, yaml_file, sort_keys=False)

    def load_config(self, config_file: str):
        '''Method to load a config file and create all the objects listed in 
        the config file

        Parameters:            
            config_file (str): Config file to load

        Returns:
            None
        '''

        if os.path.splitext(config_file)[1] not in SUPPORTED_YAML_EXTS:
            raise ValueError(
                f'The file you are loading does not seem to be a yaml file'
                f'(file must end with {" ,".join(SUPPORTED_YAML_EXTS)})')

        with open(config_file) as yaml_file:
            data_dict = yaml.safe_load(yaml_file)
            print(data_dict)
            # TODO: create system, lammps, job objects with this dict
