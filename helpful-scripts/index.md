# Index of Helpful Scripts

| Name     | Description    |
|----------|----------------|
| NEB_barrier_evolution.ipynb   | Analysis of NEB convergence (energy profile changes, errors, plot of images).     |
| neb-visualisation.ipynb       | Analysis of NEB convergence (energy profile changes, errors, plot of images), a little harder to debug.        |
| charge-neb-images.py          | Extract the images from a converged NEB calculation, apply charge, and generate Environ input file with countercharge planes.     |
| relax-structures.ipynb        | Analysis of BFGS steps of a relax calculation (does not need to be converged). Will produce a plot of the energy over BFGS steps, as well as the structures at each step.        |
| charge-relax-structure.py     | Apply a series of charges to a relaxed structure, generate relax files of those and input them into folders, and generate Environ input files with countercharge planes.     |
| run_charges.sh                | Bash file to run all the relax input files from charge-relax-structure.py in their directories.  |
