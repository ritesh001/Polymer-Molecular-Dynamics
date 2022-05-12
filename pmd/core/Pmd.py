import os
import yaml
from typing import Any, Dict, List, Union

from pmd.core.ForceField import ForceField
from pmd.core.Job import Job
from pmd.core.Lammps import Lammps
from pmd.core.Procedure import Procedure
from pmd.core.System import System
from pmd.util import Util

PRIMITIVE_TYPE_LIST = (str, int, float, bool)


# Custom yaml config file dictionary constructor
def to_yaml_dict(
        cls: Union[System, ForceField, Lammps, Procedure, Job]) -> Dict:
    return {
        # strip off the front underscore
        k.lstrip('_'): custom_class_yaml_dumper(v)
        for k, v in cls.__dict__.items()
        # only add to dict if value is not None nor System/Lammps instances
        if v is not None and not issubclass(type(v), (System, Lammps))
    }


# Custom method to dump values of non-primitive type to the yaml config file
def custom_class_yaml_dumper(v: Any) -> Any:
    # If value is a list, recursively go through each item in the list
    # Specifically, this is for the Lammps procedure list
    if isinstance(v, list):
        return [custom_class_yaml_dumper(i) for i in v]
    # If value is a non-primitive type, expand it to a dict
    elif not isinstance(v, PRIMITIVE_TYPE_LIST):
        return {str(v): to_yaml_dict(v)}
    # If value is a primitive type, simply dump it
    else:
        return v


class Pmd:
    '''Template object to perform tasks for Systems, Lammps, and Jobs 
       altogether (e.g. create files and metadata)

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
    def create(self, output_dir: str, save_config: bool = False) -> None:
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
            self.save_config(output_dir)

    @Util.build_dir
    def save_config(self, output_dir: str, config_fname: str = 'config.yaml'):
        with open(os.path.join(output_dir, config_fname), 'w') as yaml_file:
            config = {
                str(config_item): to_yaml_dict(config_item)
                for config_item in self._config_items
            }
            yaml.dump(config, yaml_file, sort_keys=False)

    def load_config(self, config_fname: str):
        with open(config_fname) as yaml_file:
            data_dict = yaml.load(yaml_file)
            print(data_dict)
            # TODO: create system, lammps, job objects with this dict
