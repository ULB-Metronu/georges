# Introduction

## Usage
No attempt is made to support python versions earlier than CPython 3.5. Jython and alternative implementation have not been tested.

### Physics module ###
`georges.physics` is a simple standalone helper module providing various relativistic conversion functions as well as range computations for protons in water.

It is automatically importe when the main `georges` module is imported.

Pull requests are encouraged, as the module growths it will be necessary to split it in autonomous pieces.

```python
import georges
georges.physics.energy_to_beta(0.1)
georges.physics.range_to_energy(4.0)
```

The georges.physics.kinematics function allows rapid and easy conversion between the different kinematics quantities used in beam physics (kinetic energy, momentum, magnetic rigidity, relativistic gamma and beta) and to proton range in water (using the IEC60601 range in water conversion).

```python
georges.physics.kinematics(energy=100)
{
 'beta': 0.43244159571922342,
 'brho': 1.5010379999999999,
 'energy': 102.33086821870768,
 'gamma': 1.109063106808982,
 'momentum': 450,
 'range': 8.0501247517130885
 }
```

```python
georges.physics.kinematics(gamma=1.2)
{
 'beta': 0.5527707983925666,
 'brho': 2.0760332515185556,
 'energy': 187.65441625999995,
 'gamma': 1.2,
 'momentum': 622.3792889875873,
 'range': 23.248200723538545
}
```

```python
georges.physics.kinematics(beta=0.4)
{
 'beta': 0.4,
 'brho': 1.365929596629474,
 'energy': 85.466688943097552,
 'gamma': 1.0910894511799618,
 'momentum': 409.49550809723888,
 'range': 5.8362733414685106
 }
```

```python
georges.physics.kinematics(range=32)
{
 'beta': 0.59219057901310057,
 'brho': 2.3000819451857275,
 'energy': 226.12911179644985,
 'gamma': 1.2410059046872013,
 'momentum': 689.54741674333184,
 'range': 32
 }
```