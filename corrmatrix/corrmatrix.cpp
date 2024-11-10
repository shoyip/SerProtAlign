#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <omp.h>

using namespace std;

// Given an alignment, the number of positions, the number
// of sequences and the position indices i and j, yields
// the value of the correlation matrix Ctilde_ij

// Read through the specified column of an alignment
// At each residue, store the value in a NxNxQxQ array

int get_aa_idx(char residue) {
    int aa_idx;
    switch (residue) {
        case '-':
            aa_idx = 0;
            break;
        case 'A':
            aa_idx = 1;
            break;
        case 'C':
            aa_idx = 2;
            break;
        case 'D':
            aa_idx = 3;
            break;
        case 'E':
            aa_idx = 4;
            break;
        case 'F':
            aa_idx = 5;
            break;
        case 'G':
            aa_idx = 6;
            break;
        case 'H':
            aa_idx = 7;
            break;
        case 'I':
            aa_idx = 8;
            break;
        case 'K':
            aa_idx = 9;
            break;
        case 'L':
            aa_idx = 10;
            break;
        case 'M':
            aa_idx = 11;
            break;
        case 'N':
            aa_idx = 12;
            break;
        case 'P':
            aa_idx = 13;
            break;
        case 'Q':
            aa_idx = 14;
            break;
        case 'R':
            aa_idx = 15;
            break;
        case 'S':
            aa_idx = 16;
            break;
        case 'T':
            aa_idx = 17;
            break;
        case 'V':
            aa_idx = 18;
            break;
        case 'W':
            aa_idx = 19;
            break;
        case 'Y':
            aa_idx = 20;
            break;
    }
    return aa_idx;
}

float get_ctilde(string filename, int i, int j, int Q) {
    int onept_count_i[Q];
    int onept_count_j[Q];
    int twopt_counts[Q][Q];
    memset(onept_count_i, 0, sizeof(onept_count_i));
    memset(onept_count_j, 0, sizeof(onept_count_j));
    memset(twopt_counts, 0, sizeof(twopt_counts));

    ifstream file(filename);

    if(!file.is_open()) {
        cerr << "Error opening file!" << endl;
    }

    string line;
    int i_aa_idx, j_aa_idx;
    int n_seqs = 0;
    while (getline(file, line)) {
        n_seqs++;
        for (size_t k=0; k<line.size(); ++k) {
            if (k==i) {
                i_aa_idx = get_aa_idx(line[i]);
                onept_count_i[i_aa_idx] = 1;
            } else if (k==j) {
                j_aa_idx = get_aa_idx(line[j]);
                onept_count_j[j_aa_idx] = 1;
            } else {
                continue;
            }
        }
        twopt_counts[i_aa_idx][j_aa_idx]++;
    }
    file.close();

    float ctilde_arg = 0.;
    float ctilde = 0.;

    for (int k=0; k<Q; ++k) {
        for (int l=0; l<Q; ++l) {
            ctilde_arg = 1. * twopt_counts[k][l] / n_seqs - 1. * onept_count_i[k] * onept_count_j[l] / (n_seqs * n_seqs);
            ctilde += ctilde_arg * ctilde_arg;
        }
    }

    return ctilde;
}

int main() {
    float ctilde_arr[693][693];
    for (int i=0; i<693; ++i) {
        for (int j=0; j<693; ++j) {
            cout << "i" << i << "j" << j << endl;
            float ctilde_val = get_ctilde("aln_short_short.txt", i, j, 21);
            ctilde_arr[i][j] = ctilde_val;
        }
    }

    ofstream ctilde_file;
    ctilde_file.open("ctilde.csv");
    for (int i=0; i<693; ++i) {
        for (int j=0; j<693; ++j) {
            ctilde_file << setprecision(6) << ctilde_arr[i][j] << ",";
        }
        ctilde_file << endl;
    }
    ctilde_file.close();
    return 0;
}
