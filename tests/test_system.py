import pytest

from pmd.core import EMC, PSP, SolventSystem, System


@pytest.fixture
def test_data():
    test_builder = PSP('opls-lbcc')

    syst = System('*CC*',
                  builder=test_builder,
                  density=0.5,
                  natoms_total=500,
                  natoms_per_chain=100)

    solventsyst = SolventSystem(smiles='*CC*',
                                solvent_smiles='C1CCCCC1',
                                ru_nsolvent_ratio=0.098469007,
                                builder=test_builder,
                                density=0.8,
                                natoms_total=5000,
                                natoms_per_chain=150)
    return {
        'builder': test_builder,
        'system': syst,
        'solventsystem': solventsyst,
    }


def test_system_initialization(test_data):
    syst = test_data['system']
    assert syst.smiles == '*CC*'
    assert syst.builder == test_data['builder']


def test_system_update(test_data):
    syst = test_data['system']
    new_builder = EMC('pcff')
    syst.smiles = '*CC(*)CC'
    syst.builder = new_builder
    assert syst.smiles == '*CC(*)CC'
    assert syst.builder == new_builder


def test_solventsystem_initialization(test_data):
    solv_syst = test_data['solventsystem']
    assert solv_syst.smiles == '*CC*'
    assert solv_syst.builder == test_data['builder']


def test_solventsystem_update(test_data):
    solv_syst = test_data['solventsystem']
    solv_syst.smiles = '*CC(*)CC'
    assert solv_syst.smiles == '*CC(*)CC'


def test_system_write_data(tmp_path, test_data):
    # PSP not installed
    syst = test_data['system']
    with pytest.raises(ImportError):
        syst.write_data()

    d = tmp_path / "result"
    syst.builder = EMC('pcff')
    syst.write_data(d)
    assert len(list(tmp_path.iterdir())) == 1
