from pmd.postprocessing.TrajectoryReader import read_lammpstrj, read_lammpstrj_by_type


def test_read_lammpstrj(data_path):
    f = data_path / 'equil.lammpstrj'
    r, ir, timestep, box_bounds, id2type, id2mol, mol2ids = read_lammpstrj(f)

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


def test_read_lammpstrj_by_type(data_path):
    f = data_path / 'equil.lammpstrj'
    (r, ir, timestep, box_bounds, id2type,
     id2nidex) = read_lammpstrj_by_type(f, types=[1, 2])

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