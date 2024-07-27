#include "molecule.h"
#include "masses.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
using namespace std;


int main(int argc, char* argv[])
{
    if (argc != 2){
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    Molecule *mol = nullptr;

    try {
        mol = new Molecule(argv[1], 0);
        cout << "input geometry:"  << endl;
        mol->print_geom();
    } catch (const exception &e) {
        cerr << e.what() << endl;
        delete mol;
        return 1;
    }

    cout << "Interatomic distances (bohr):\n";
    for(int i=0; i < mol->natom; i++)
        for(int j=0; j < i; j++)
            printf("%d %d %8.5f\n", i, j, mol->bond(i,j));


    // print angles
    cout << "\nBond angles \n";
    for (int i = 0; i < mol->natom; ++i){
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                if (mol->bond(i,j) < 4.0 && mol->bond(j,k) < 4.0)
                    printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol->angle(i, j, k)*(180.0/acos(-1.0)));
            }
        }
    }

    // finding the center of mass
    double M = 0.0;
    for (int i = 0; i < mol->natom; i++)
        M += masses.at(mol->zvals[i]);

    double xcm=0.0;
    double ycm=0.0;
    double zcm=0.0;
    double mi;
    for(int i=0; i < mol->natom; i++) {
        mi = masses.at((int) mol->zvals[i]);
        xcm += mi * mol->geom[i][0];
        ycm += mi * mol->geom[i][1];
        zcm += mi * mol->geom[i][2];
    }
    xcm /= M;
    ycm /= M;
    zcm /= M;

    printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", xcm, ycm, zcm);
    // translate molecule to center of mass
    mol->translate(-xcm, -ycm, -zcm);

    // define tensor of inertia
    Matrix I(3, 3);
    for(int i=0; i < mol->natom; i++) {
        mi = masses.at((int) mol->zvals[i]);
        I(0,0) += mi * (mol->geom[i][1]*mol->geom[i][1] + mol->geom[i][2]*mol->geom[i][2]);
        I(1,1) += mi * (mol->geom[i][0]*mol->geom[i][0] + mol->geom[i][2]*mol->geom[i][2]);
        I(2,2) += mi * (mol->geom[i][0]*mol->geom[i][0] + mol->geom[i][1]*mol->geom[i][1]);
        I(0,1) -= mi * mol->geom[i][0]*mol->geom[i][1];
        I(0,2) -= mi * mol->geom[i][0]*mol->geom[i][2];
        I(1,2) -= mi * mol->geom[i][1]*mol->geom[i][2];
    }
    I(1,0) = I(0, 1);
    I(2,0) = I(0, 2);
    I(2,1) = I(1, 2);

    cout << "\n Moment of inertia tensor in [amu*bohr^2]:\n";
    cout << I << endl;

    // find principal moments, our matrix is symmetric => self adjoint
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Vector evals = solver.eigenvalues();  
    cout << "\nPrincipal moments of inertia" << endl;
    cout << evals << endl;

    // classify rigid rotors https://en.citizendium.org/wiki/Classification_of_rigid_rotors
    if(mol->natom == 2) 
        cout << "Molecule is diatomic" << endl;
    else if (evals[0] < 1e-4) 
        cout << "Molecule is linear" << endl;
    else if (fabs(evals[0] - evals[1]) < 1e-4 && fabs(evals[1] - evals[2]) < 1e-4)
        cout << "Molecule is a spherical top" << endl;
    else if (fabs(evals[0] - evals[1]) < 1e-4 && fabs(evals[1] - evals[2]) > 1e-4)
        cout << "Molecule is an oblate symmetric top" << endl;
    else if (fabs(evals[0] - evals[1]) > 1e-4 && fabs(evals[1] - evals[2]) < 1e-4)
        cout << "Molecule is a prolate symmetric top" << endl;
    else cout << "Molecule is an asymmetric top" << endl;

    
    // compute the rotational constants 
    double pi = acos(-1.0);
    double conv = 0;
    conv = 6.6260755e-34/(8.0 * pi * pi);
    conv /= 1.6605402e-27 * 0.529177249e-10 * 0.529177249e-10;
    conv *= 1e-6;
    cout << "\nRotational constants in [MHz]:" << endl;
    cout << "\tA = " << conv/evals[0] << "\t B = " << conv/evals[1] << "\t C = " << conv/evals[2] << endl;


    return 0;
}