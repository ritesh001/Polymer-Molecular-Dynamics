import pytest

from pmd.core import Torque


@pytest.mark.parametrize(
    "job, job_fname",
    [
        (Torque(run_lammps='lmp.in',
                jobname='PE_equilibration',
                project='GT-rramprasad3-CODA20',
                nodes=1,
                ppn=24,
                walltime='24:00:00'), 'job.pbs'),
    ],
)
def test_job_write(data_path, tmp_path, job, job_fname):
    d = tmp_path / "result"
    expected_output = data_path / job_fname
    actual_output = d / job_fname
    job.write_job(d)

    assert actual_output.read_text() == expected_output.read_text()
