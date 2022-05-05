import os
import json
from typing import List, Union

from pmd.core.Job import Job
from pmd.core.Lammps import Lammps
from pmd.core.System import System


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

    def create(self, output_dir: str, metadata: bool = False) -> None:
        if self._system:
            self._system.write_data(output_dir)
        if self._lammps:
            for lmp in self._lammps:
                lmp.write_lammps(output_dir)
        if self._job:
            for job in self._job:
                job.write_job(output_dir)

        if metadata:
            self.create_metadata(output_dir)

    def create_metadata(self, output_dir: str, metadata_fname='metadata.json'):
        with open(os.path.join(output_dir, metadata_fname), 'w') as json_file:
            metadata = {
                'system': str(self._system),
                'lammps': str(self._lammps),
                'job': str(self._job),
            }
            json.dump(metadata, json_file)

    def load_metadata(self, metadata_fname: str):
        with open(metadata_fname) as json_file:
            data_dict = json.loads(json_file)
            print(data_dict)
            # TODO: create system, lammps, job objects with this dict
