# Density Functional Theory
This is the repository for the master work at the University of Oslo.

## Structure
- Code
This folder contains all of the code used for pre and post analysis of data
concerning DFT calculations.

  - VASP
    - .. (Holds all the VASP input files that were used in the project)
  - olem_scripts
    - .. (Contains scripts that can be used in extracting data from VASP
      output files)
  - pymatgen
    - .. (A python package that was used to plot the data from VASP output files)
  - python_scripts
    - .. (Contains scripts that were used in this project to construct and
      extract data)

- Data
This folder contains files for preparation of data, such as atom positions for DFT
calculations.

  - Layer
    - .. (All POSCAR files that were constructed regarding the 2D layered structures)

- Results
This folder contains all the results produced from the post analysis from the
DFT calculations.

  - Bulk
    - .. (All output files from VASP that concerns the bulk structures)
  - Figures
    - .. (All PDF's from the notebooks)
  - Layers
    - .. (All output files from VASP that concerns the 2D layered structures)
  - Notebooks
    - .. (This is where the jupyter notebooks are, where all the reusults are
      displayed)
  - atom_energy
    - .. (Energy calculations from VASP that concerns the individual atoms of
      the constituent atoms)
  - solid_energy
    - .. (Energy calculations from VASP that concerns the solid structures of
      the constituent atoms)

## How to run the files
To run the notebooks and the python scripts, one must utilize the virtual
environment `pipenv`, therefore one must first install the package `pipenv`. 
After installing $pipenv$, one must initialize the virtual environment, and this can be done by simply inputting:

```
$ pipenv shell 
$ pipenv sync
```
 

## Resources:

### Links
- https://pymatgen.org/
  - https://pymatgen.org/modules.html
    - https://pymatgen.org/pymatgen.electronic_structure.html
      - https://pymatgen.org/pymatgen.electronic_structure.plotter.html
    - https://pymatgen.org/pymatgen.io.vasp.html
       - https://pymatgen.org/pymatgen.io.vasp.outputs.html

- vaspWiki
  - https://www.vasp.at/wiki/index.php/Category:Output_Files
