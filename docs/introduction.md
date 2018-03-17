# Introduction

## Usage
No attempt is made to support python versions earlier than CPython 3.5. Jython and alternative implementation have not been tested.

### Physics module ###
`georges.physics` is a simple standalone helper module providing various relativistic conversion functions as well as range computations for protons in water.

It is automatically importe when the main `georges` module is imported.

Pull requests are encouraged, as the module growths it will be necessary to split it in autonomous pieces.

```
import georges
georges.physics.energy_to_beta(0.1)
georges.physics.range_to_energy(4.0)
```
