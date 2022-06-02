import pytest

from pmd.core import GAFF2, OPLS, SolventSystem, System

# TODO: put into a pytest.mark.parametrize
test_force_field = GAFF2()
syst = System('*CC*',
              force_field=test_force_field,
              density=0.8,
              natoms_total=5000,
              natoms_per_chain=150)
solv_syst = SolventSystem(smiles='*CC*',
                          solvent_smiles='C1CCCCC1',
                          ru_nsolvent_ratio=0.098469007,
                          force_field=test_force_field,
                          density=0.8,
                          natoms_total=5000,
                          natoms_per_chain=150)


def test_system_initialization():
    assert syst.smiles == '[*]CC[*]'
    assert syst.force_field == test_force_field


def test_system_update():
    new_force_field = OPLS()
    syst.smiles = '*CC(*)CC'
    syst.force_field = new_force_field
    assert syst.smiles == '[*]CC([*])CC'
    assert syst.force_field == new_force_field


def test_solventsystem_initialization():
    assert solv_syst.smiles == '[*]CC[*]'
    assert solv_syst.force_field == test_force_field


def test_solventsystem_update():
    new_force_field = OPLS()
    solv_syst.smiles = '*CC(*)CC'
    solv_syst.force_field = new_force_field
    assert solv_syst.smiles == '[*]CC([*])CC'
    assert solv_syst.force_field == new_force_field


def test_system_write_data():
    # TODO: properly test write_data
    with pytest.raises(ImportError):
        syst.write_data()
