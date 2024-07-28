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
    if (argc != 3){
    cerr << "Usage: " << argv[0] << " <molecule.dat> <hessian>" << endl;
    return 1;
    }

    // create molecule
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

    // read hessian
    Matrix H(3*mol->natom, 3*mol->natom);
    string **indices = new string* [3*mol->natom];
    for (int i=0; i<3*mol->natom; ++i)
        indices[i] = new string[3*mol->natom];
    // open fstream
    ifstream hesfile(argv[2]);
    int nat;
    // read number of columns in compressed hessian
    hesfile >> nat;
    printf("Hessian file has %d atoms\n", nat);
    // read file
    for (int i=0; i<3*mol->natom; i++)
        for (int j=0; j<mol->natom; ++j) {
            hesfile >> H(i, 3*j);
            hesfile >> H(i, 3*j+1);
            hesfile >> H(i, 3*j+2);
            indices[i][3*j] = to_string(i / 3) + "," + to_string(3*j / 3);
            indices[i][3*j+1] = to_string(i / 3) + "," + to_string(3*j / 3);
            indices[i][3*j+2] = to_string(i / 3) + "," + to_string(3*j / 3);
        }
    hesfile.close();

    // https://www.researchgate.net/figure/Structure-of-the-Hessian-matrix-H-for-a-system-containing-N-atom-4-atoms-The-atomic_fig1_358914864
    cout << "\nindiced of hessian\n";
    for (int i=0; i<3*mol->natom; i++) {
        for (int j=0; j<3*mol->natom; j++) {
            cout << indices[i][j] << " ";
        }
        cout << endl;
    }
    cout << "\nReaded Hessian:\n" << H << std::endl;

    // mass-weighten the hessian
    for (int i=0; i<3*mol->natom; ++i)
        for (int j=0; j<3*mol->natom; ++j) {
            H(i, j) /= sqrt(masses.at(mol->zvals[i/3]) * masses.at(mol->zvals[j/3]));
        }
    cout << "\nMass-weighten Hessian:\n" << H << std::endl;

    // compute eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H);
    Matrix eigenvec = solver.eigenvectors();
    Vector eigenval = solver.eigenvalues();
    cout << "\nEigenvalues of Hessian:\n";
    cout << eigenval << endl;

    // convert to cm-1
    double conv = 5140.489580228;
    for (int i=0; i<eigenval.rows(); ++i)
        eigenval[i] = conv * sqrt(eigenval[i]);
    cout << "\nVibrational frequencies in cm-1:\n" << eigenval;

    for (int i=0; i<3*mol->natom; ++i)
        delete[] indices[i];
    delete[] indices;
    return 0;
}