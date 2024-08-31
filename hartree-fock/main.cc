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

int get_nrows_2d(const string &filename) {
    ifstream input(filename);
    if (!input.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        throw runtime_error("Error reading 1e integral file");
    }
    int n = 0;
    int nn;
    string line;

    while (input >> nn >> line >> line) {
        if (nn>n) {
            n = nn;
        }
    }

    input.close();
    return n;
}

int read_2d(const string &filename, Matrix &M, int n) {
    // calculate number of elements to read
    int N;
    int i, j;
    N = n * (n + 1) / 2;
    ifstream input(filename);
    if (!input.is_open()) {
        cerr << "Error opening " << filename << endl;
    }

    for (int k=0; k<N; ++k) {
        input >> i >> j;
        input >> M(i-1, j-1);
        //cout << M(i-1, j-1);
        if (i != j) {
            M(j-1, i-1) = M(i-1, j-1);
        }
    }
    input.close();
    return 0;
}

int read_4d(const string &filename, double**** V) {
    ifstream input(filename);
    if (!input.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        throw runtime_error("Error reading 1e integral file");
    }

    int i, j, k, l;
    double val;

    while (input >> i >> j >> k >> l >> val) {
        --i;
        --j;
        --k;
        --l;
        V[i][j][k][l] = val; // 1
        // restore by symmetry
        V[j][i][k][l] = val; // 2
        V[i][j][l][k] = val; // 3
        V[j][i][l][k] = val; // 4
        V[k][l][i][j] = val; // 5   
        V[l][k][i][j] = val; // 6
        V[k][l][j][i] = val; // 7
        V[l][k][j][i] = val; // 8
    }

    return 0;
}



int main(int argc, char* argv[])
{
    double conv_e = 1e-10;
    double conv_d = 1e-9;
    int max_iter = 100;
    
    if (argc != 2){
        cerr << "Usage: " << argv[0] << " <./example/folder/>" << endl;
        return 1;
    }
    string folder;
    folder = argv[1];
    //cout << folder << endl;

    // create molecule
    Molecule *mol = nullptr;

    try {
        mol = new Molecule((folder + "geom.dat").c_str(), 0);
        cout << "input geometry:"  << endl;
        mol->print_geom();
    } catch (const exception &e) {
        cerr << e.what() << endl;
        delete mol;
        return 1;
    }

    // read enuc
    ifstream input(folder + "enuc.dat");
    if (!input.is_open()) {
        std::cerr << "Error opening file: " << folder + "enuc.dat" << std::endl;
        return 1;
    }
    double Enuc;
    input >> Enuc;
    input.close();

    // get number of rows in upper triangular form
    int n;
    try {
        n = get_nrows_2d(folder + "s.dat");
    }
    catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return 1;
    }

    // read overlap integrals
    Matrix S(n, n);
    read_2d(folder + "s.dat", S, n);
    // read kinetic integrals
    Matrix T(n, n);
    read_2d(folder + "t.dat", T, n);
    // read pot integrals
    Matrix V(n, n);
    read_2d(folder + "v.dat", V, n);

    Matrix H_core(n, n);
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            H_core(i, j) = V(i, j) + T(i, j);
        }
    }

    // compute indices for 2e integrals
    // Allocate memory for 4D array
    double ****int2e = new double***[n];
    for (int i=0; i<n; ++i) {
        int2e[i] = new double**[n];
        for (int j=0; j<n; ++j) {
            int2e[i][j] = new double*[n];
            for (int k=0; k<n; ++k) {
                int2e[i][j][k] = new double[n];
                for (int l=0; l<n; ++l) {
                    int2e[i][j][k][l] = 0.00000000;
                }
            }
        }
    } 
    // read integrals from file
    read_4d(folder + "eri.dat", int2e);

    // Calculate number of occupied orbitals
    int nelec = 0;
    for (int i=0; i<mol->natom; ++i) {
        nelec += mol->zvals[i];
    }
    int nocc = nelec / 2;

    // diagonalize overlap matrix
    Eigen::SelfAdjointEigenSolver<Matrix> solver(S); // use SelfAdjoint property to avoid complex eigenvalues = types
    Vector D = solver.eigenvalues();
    Matrix L = solver.eigenvectors();

    // Compute the inverse square root of the eigenvalues
    Vector D_inv_sq = D.array().inverse().sqrt();

    // Convert D_inv_sq to a diagonal matrix
    Matrix D_inv_sq_diag = D_inv_sq.asDiagonal();
    Matrix S_inv_sq = L * D_inv_sq_diag * L.transpose();

    //cout << S_inv_sq << endl;

    // Build initial Fock matrix in orthonormal AO basis with H_core as guess

    Matrix F_init = S_inv_sq.transpose() * H_core * S_inv_sq;

    
    // diagonalize Fock matrix to get orbitals and their energies
    solver.compute(F_init);
    Vector e_init = solver.eigenvalues();
    Matrix F_eigenvec = solver.eigenvectors();

    // Transform eigenvectors to original AO basis
    Matrix C0 = S_inv_sq * F_eigenvec;
    
    // form density matrix
    Matrix Dens(n, n);
    Dens.setZero();
    for (int m=0; m<nocc; ++m) {
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                Dens(i, j) += C0(i, m) * C0(j, m);
            }
        }
    }
    
    // Compute SCF energy
    double Eelec = 0.0;

    // Iterate over all elements of the matrices
    for (int mu = 0; mu < Dens.rows(); ++mu) {
        for (int nu = 0; nu < Dens.cols(); ++nu) {
            Eelec += Dens(mu, nu) * (H_core(mu, nu) + F_init(mu, nu));
        }
    }
    double Etotal = Eelec + Enuc;
    printf("Iteration 0, energy: %12.6f\n", Etotal);

    Matrix Dens_old = Dens;
    double Eelec_old = Eelec;
    Matrix F = H_core;
    Matrix C(F.rows(), F.rows());
    int num_iter = 0;
    Vector orb_e(F.rows());

    while (true) {
        ++num_iter;
        // Form new Fock matrix
        F = H_core;
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                for (int k=0; k<n; ++k) {
                    for (int l=0; l<n; ++l) {
                        F(i, j) += Dens(k, l) * (2 * int2e[i][j][k][l] - int2e[i][k][j][l]);
                    }
                }
            }
        }

        // Form new density matrix
        F = S_inv_sq.transpose() * F * S_inv_sq;
        solver.compute(F);
        orb_e = solver.eigenvalues();  // diagonalize Fock matrix in orthogonal basis
        Matrix C_prime = solver.eigenvectors();
        C = S_inv_sq * C_prime;
        Dens.setZero();
        for (int m=0; m<nocc; ++m) {
            for (int i=0; i<n; ++i) {
                for (int j=0; j<n; ++j) {
                    Dens(i, j) += C(i, m) * C(j, m);
                }
            }
        }
        // Compute new SCF energy
        Eelec = 0.0;
        for (int mu = 0; mu < Dens.rows(); ++mu) {
            for (int nu = 0; nu < Dens.rows(); ++nu) {
                Eelec += Dens(mu, nu) * (H_core(mu, nu) + F(mu, nu));
            }
        }
        
        Etotal = Eelec + Enuc;

        // Check for convergence
        double rmsd = 0;
        for (int mu = 0; mu < Dens.cols(); ++mu) {
            for (int nu = 0; nu < Dens.rows(); ++nu) {
                rmsd += pow(Dens(mu, nu) - Dens_old(mu, nu), 2);
            }
        }
        rmsd = sqrt(rmsd);
        printf("Iteration %d, energy: %16.10f, density matrix rmsd: %16.10f\n", num_iter, Etotal, rmsd);

        if (abs(Eelec_old-Eelec) < conv_e && rmsd < conv_d) {
            break;
        }
        Dens_old = Dens;
        Eelec_old = Eelec;

        // Exit if not converged        
        if (num_iter == max_iter) {
            cout << "Self consisted field did not converged" << endl;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        delete[] int2e[i][j][k];
                    }
                    delete[] int2e[i][j];
                }
                delete[] int2e[i];
            }
            delete[] int2e;
            return 1;
        }
    }

    cout << "Succesfull convergence in " << num_iter << " iterations" << endl;
    cout << orb_e << endl;
    
    // transform Fock matrix to MO basis
    Matrix F_mo(F.rows(), F.rows());
    F_mo.setZero();
    // for (int i = 0; i<F_mo.rows(); ++i) {
    //     for (int j = 0; j < F_mo.rows(); ++j) {
    //         for (int mu = 0; mu < F_mo.rows(); ++mu) {
    //             for (int nu = 0; nu < F_mo.rows(); ++nu) {
    //                 F_mo(i, j) += C(j, mu) * C(i, nu) * F(mu, nu);
    //             }
    //         }
    //     }
    // }
    F_mo = C.transpose() * F * C;
    cout << F << endl;
















    // deallocate 2 electron integrals
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                delete[] int2e[i][j][k];
            }
            delete[] int2e[i][j];
        }
        delete[] int2e[i];
    }
    delete[] int2e;
    return 0;
}



