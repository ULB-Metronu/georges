import pytest

import pymadx

@pytest.fixture()
def reference_beam():
    return pymadx.Beam.Beam(particletype='e-',
                            energy=1.0,
                            distrtype='reference')

def test_SetEnergy(reference_beam):
    reference_beam.SetEnergy(10.0)
    assert reference_beam.GetItemStr('energy') == '10.0'

@pytest.mark.parametrize('test_input, expected',
                         [('e-', 'electron'),
                          ('e+', 'positron'),
                          ('proton', 'proton')])
def test_SetParticleType(reference_beam, test_input, expected):
    reference_beam.SetParticleType(test_input)
    assert reference_beam['particle'] == expected

def test_SetParticleType_raises_exception(reference_beam):
    with pytest.raises(ValueError):
        reference_beam.SetParticleType('nonsense')
