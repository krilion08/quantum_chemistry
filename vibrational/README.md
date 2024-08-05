## compile:
>g++ -I. -I./Eigen -o read.x read_molecule.cc molecule.cc masses.cc -std=c++11

## test

>./vibrations.x examples/<molecule>.dat examples/<molecule_hessian>.txt