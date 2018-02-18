# Fermi-Eyges

The `Fermi-Eyges` module is a reimplementation of the Fermi-Eyges transport framework, largely based on the work and publications of Bernard Gottschalk.

The module is composed of:
 - *Materials database* (`materials_db.py`): reads, loads and provides an interface to data (stopping power, radiation lengths, etc.) for a large set of materials commonly used in protontherapy. When using the Fermi-Eyges module, a database must be instanciated and passed around to the other components. The database reads stopping power and range tables from the p-star database format;
 - *Stopping power* (`stopping.py`): computes the various range-energy direct and inverse relationships from measured range stopping power data tables;
 - *MCS* (multiple Coulomb scattering) (`mcs.py`): 
 

# Typical use case


```
import georges.fermi
db = georges.fermi.MaterialsDB()
georges.fermi.get_range_from_energy(material='water',
                                   energy=70,
                                   db=db)
                                   
georges.fermi.get_energy_from_range(material='graphite',
                                   r=10,
                                   db=db
                                   )
                                   
georges.fermi.compute_fermi_eyges(material='graphite', 
                                  energy=228.15, 
                                  thickness=18.5, 
                                  db=db, 
                                  T=georges.fermi.DifferentialMoliere)
                                  
georges.fermi.residual_energy(material='graphite', thickness=1, k_in=50, db=db)                                                                                                       
```
 