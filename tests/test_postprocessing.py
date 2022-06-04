import pytest

from pmd.entry import analyze
from pmd.postprocessing.Analysis import calculate_MSD, calculate_Tg
from pmd.postprocessing.TrajectoryReader import (read_lammpstrj,
                                                 read_lammpstrj_by_type)


@pytest.fixture
def test_data(data_path):
    return {
        'trajectory_file': data_path / 'equil.lammpstrj',
        'Tg_result_file': data_path / 'temp_vs_density.txt'
    }


def test_read_lammpstrj(test_data):
    (r, ir, timestep, box_bounds, id2type, id2mol,
     mol2ids) = read_lammpstrj(test_data['trajectory_file'])

    assert len(r) == 2  # 2 timestep
    assert len(r[0][1:]) == 1248  # 1248 atoms
    assert len(r[0][0]) == 3  # x, y, z

    assert len(ir) == 2  # 2 timestep
    assert len(ir[0][1:]) == 1248  # 1248 atoms
    assert len(ir[0][0]) == 3  # x, y, z

    assert len(timestep) == 2  # 2 timestep

    assert len(box_bounds) == 2  # 2 timestep
    assert len(box_bounds[0]) == 3  # x, y, z
    assert len(box_bounds[0][0]) == 2  # lower, upper

    assert len(id2type[1:]) == 1248  # 1248 atoms
    assert len(id2mol[1:]) == 1248  # 1248 atoms
    assert len(mol2ids[1:]) == 12  # 12 polymer chains


def test_read_lammpstrj_skip_beginning(test_data):
    (r, ir, timestep, box_bounds, id2type, id2mol,
     mol2ids) = read_lammpstrj(test_data['trajectory_file'], skip_beginning=1)

    assert len(r) == 1  # 1 timestep
    assert len(r[0][1:]) == 1248  # 1248 atoms
    assert len(r[0][0]) == 3  # x, y, z

    assert len(ir) == 1  # 1 timestep
    assert len(ir[0][1:]) == 1248  # 1248 atoms
    assert len(ir[0][0]) == 3  # x, y, z

    assert len(timestep) == 1  # 1 timestep

    assert len(box_bounds) == 1  # 1 timestep
    assert len(box_bounds[0]) == 3  # x, y, z
    assert len(box_bounds[0][0]) == 2  # lower, upper

    assert len(id2type[1:]) == 1248  # 1248 atoms
    assert len(id2mol[1:]) == 1248  # 1248 atoms
    assert len(mol2ids[1:]) == 12  # 12 polymer chains


def test_read_lammpstrj_by_type(test_data):
    (r, ir, timestep, box_bounds, id2type,
     id2nidex) = read_lammpstrj_by_type(test_data['trajectory_file'],
                                        types=[1, 2])

    assert len(r) == 2  # 2 timestep
    assert len(r[0][1:]) == 1248  # 1248 atoms
    assert len(r[0][0]) == 3  # x, y, z

    assert len(ir) == 2  # 2 timestep
    assert len(ir[0][1:]) == 1248  # 1248 atoms
    assert len(ir[0][0]) == 3  # x, y, z

    assert len(timestep) == 2  # 2 timestep

    assert len(box_bounds) == 2  # 2 timestep
    assert len(box_bounds[0]) == 3  # x, y, z
    assert len(box_bounds[0][0]) == 2  # lower, upper

    assert len(id2type[1:]) == 1248  # 1248 atoms
    assert len(id2nidex[1:]) == 1248  # 1248 atoms


def test_read_lammpstrj_by_type_skip_beginning(test_data):
    (r, ir, timestep, box_bounds, id2type,
     id2nidex) = read_lammpstrj_by_type(test_data['trajectory_file'],
                                        types=[1, 2],
                                        skip_beginning=1)

    assert len(r) == 1  # 1 timestep
    assert len(r[0][1:]) == 1248  # 1248 atoms
    assert len(r[0][0]) == 3  # x, y, z

    assert len(ir) == 1  # 1 timestep
    assert len(ir[0][1:]) == 1248  # 1248 atoms
    assert len(ir[0][0]) == 3  # x, y, z

    assert len(timestep) == 1  # 1 timestep

    assert len(box_bounds) == 1  # 1 timestep
    assert len(box_bounds[0]) == 3  # x, y, z
    assert len(box_bounds[0][0]) == 2  # lower, upper

    assert len(id2type[1:]) == 1248  # 1248 atoms
    assert len(id2nidex[1:]) == 1248  # 1248 atoms


def test_calculate_msd(test_data):
    (r, ir, timestep, box_bounds, id2type, id2mol,
     mol2ids) = read_lammpstrj(test_data['trajectory_file'])
    msd_dict = calculate_MSD(r, ir, box_bounds, id2type)

    assert len(msd_dict) == 2  # 2 atom types
    assert len(msd_dict[1]) == 2  # 2 time frames


def test_calculate_Tg(test_data):
    Tg = calculate_Tg(test_data['Tg_result_file'])

    assert round(Tg, 1) == 295.2


def test_calculate_Tg_cli(caplog, test_data):
    Tg = analyze.main([str(test_data['Tg_result_file']), '-p', 'Tg'])

    assert 'Glass transition temperature' in caplog.text
    assert round(Tg, 1) == 295.2
