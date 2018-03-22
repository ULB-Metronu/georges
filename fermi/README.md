## Fermi-Eyges

The `Fermi-Eyges` module is a reimplementation of the Fermi-Eyges transport framework, largely based on the work and publications of Bernard Gottschalk.

The module is composed of:
 - *Materials database* (`materials_db.py`): reads, loads and provides an interface to data (stopping power, radiation lengths, etc.) for a large set of materials commonly used in protontherapy. When using the Fermi-Eyges module, a database must be instanciated and passed around to the other components. The database reads stopping power and range tables from the p-star database format;
 - *Stopping power* (`stopping.py`): computes the various range-energy direct and inverse relationships from measured range stopping power data tables;
 - *MCS* (multiple Coulomb scattering) (`mcs.py`): implements various scattering angle models (Gottschalk's `DifferentialMoliere` being the default);
 - *Fermi-Eyges* (`fermi_eyges.py`): computes the transport integrales (A0, A1, A2 and B) for a given material, thickness and incident energy. It also returns the residual energy to allow easy chaining through multiple slabs;
 - *Propagation* (`propagation.py`): propagation of a beam through a mutli-material beamline following the Fermi-Eyges transport theory.
 

### Typical use case


```python
import georges.fermi
db = georges.fermi.MaterialsDB()
georges.fermi.get_range_from_energy(material='water',
                                   energy=70,
                                   db=db)
                                   
georges.fermi.get_energy_from_range(material='graphite',
                                   r=10,
                                   db=db
                                   )
                                   
georges.fermi.residual_energy(material='graphite', thickness=1, k_in=50, db=db)                                                                                                       

georges.fermi.compute_fermi_eyges(material='graphite', 
                                  energy=228.15, 
                                  thickness=18.5, 
                                  db=db, 
                                  T=georges.fermi.DifferentialMoliere)
                                  
bl = georges.Beamline([
    {
        'NAME': 'shifter',
        'AT_ENTRY': 0.0,
        'TYPE': 'slab',
        'MATERIAL': 'beryllium',
        'LENGTH': 0.05,
    },
    {
        'NAME': 'scatterer',
        'AT_ENTRY': 0.5,
        'TYPE': 'slab',
        'MATERIAL': 'graphite',
        'LENGTH': 0.1,
    },
    {
        'NAME': 'scatterer2',
        'AT_ENTRY': 0.75,
        'TYPE': 'slab',
        'MATERIAL': 'water',
        'LENGTH': 0.01,
    },
    {
        'NAME': 'isocenter',
        'AT_ENTRY': 1.0,
        'TYPE': 'mp'
    }
], name="NOZZLE")

pbs = georges.fermi.propagate(
                        line=bl, 
                        db=db,
                        beam={
                            'energy': 230,
                            'A0': 0,
                            'A1': 0,
                            'A2': 0,
                        },
                        gaps='air'                                  
```

Plotting support is provided in the `georges/plotting` module for the visualization of scattering beamlines:

```python
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
georges.plotting.scattering(ax, pbs)
```
 