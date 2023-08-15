# Description
The scripts and notebook provided here provide a systematic way to directly calculate neutron reflectivity profiles from molecular dynamics simulations for
bilayer systems

# Dependences
Refnx python package

# How to use
Obtain the SLD for a bilayer using the provided script `sld_multiple_deutration_amber.py`. Please recenter the output of the script such that the electron density is symmetric around 0.

We assume you have 4 experimental reflectivity profiles for your systems (i) reflectivity profile from silica-D<sub>2</sub>O system, from a supported bilayer in (ii) D<sub>2</sub>O
   (iii) H<sub>2</sub>O (iv) Contrast-matched silicon (CMSi) contrast (~38% D<sub>2</sub>O).

The notebook uses silica reflectivity to obtain the Substrate SLD. The obtained silica substrate SLD is merged into the system SLD at a given position automatically determined
by the script. 

Using the experimental reflectivity in the three solvent contrast, the script varies $\alpha$, the water fraction at the silica-bilayer interface, and $\gamma$, the fraction of area
covered by water patches.
# Example files
The `example-data` has typical neutron scattering length density profiles for the MC3H-DOPC system (refer to the reference) and example reflectivity profiles.

# References
Please find the details at https://pubs.rsc.org/en/content/articlelanding/2023/nr/d3nr00987d and cite the reference: Nanoscale, 2023,15, 11647-11656 
