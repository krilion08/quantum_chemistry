#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>

using namespace std;

void Molecule::print_geom()
{
    for(int i=0; i < natom; i++)
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z)
{
    for (int i=0; i < natom; i++){
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

void Molecule::rotate(double phi){ }

double Molecule::bond(int atom1, int atom2)
{
    return sqrt(pow(geom[atom1][0]-geom[atom2][0], 2) + pow(geom[atom1][1]-geom[atom2][1], 2) + pow(geom[atom1][2]-geom[atom2][2], 2));
}

double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] - geom[b][cart])/bond(a, b);
}

double Molecule::angle(int a, int b, int c)
{
    return acos(unit(0,b,a) * unit(0,b,c) + unit(1,b,a) * unit(1,b,c) + unit(2,b,a) * unit(2,b,c));
}

double Molecule::torsion(int atom1, int atom2, int atom3, int atom4)
{
    return 0.0;
}

Molecule::Molecule(const char *filename, int q)
{   
    charge = q;

    std::ifstream input(filename);
    if (!input) {
        throw runtime_error("Error opening input file: " + string(filename));
    }

    input >> natom;
    if (input.fail()){
        throw runtime_error("Error reading natoms");
    }
    zvals = new int[natom];
    geom = new double* [natom];
    for (int i = 0; i < natom; i++){
        geom[i] = new double[3];
    }
    for (int i = 0; i < natom; i++) {
        input >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
        if (input.fail()){
            throw runtime_error("Error reading data for atom" + to_string(i+1));
        }
    }

    input.close();

}

Molecule::~Molecule()
{
  delete[] zvals;
  for(int i=0; i < natom; i++)
    delete[] geom[i];
  delete[] geom;
}