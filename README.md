# VASP analyzer

> :warning: This is a WIP (work in progress), no guarantee that this even works properly as of right now!

Small package to load vasprun.xml files and get the data in a format to perform analysis on the spin-splitting of antiferromagnetic materials. Features:
- Load sections of vasprun.xml that contains INCAR, EIGENVAL, DOSCAR and PROCAR information
- generate Bandstructure and DOS plots
- calculate spin splitting
