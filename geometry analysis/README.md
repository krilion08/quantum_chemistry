## compile:
go to /src
>g++ -I. -I./Eigen -o read.x read_molecule.cc molecule.cc masses.cc -std=c++11

## test

./read.x ../examples/<molecule>.dat