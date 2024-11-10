#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

using namespace std;

const int AA_LOOKUP[128] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 0-15
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 16-31
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,  // 32-47  ('-' is at 45)
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 48-63
    0,1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,  // 64-79   ('A' is at 65)
    13,14,15,16,17,0,18,19,0,20,0,0,0,0,0,0,  // 80-95
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 96-111
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   // 112-127
};

// Lookup table for amino acid indices
//const int AA_LOOKUP[128] = {
//    ['A'] = 1,  ['C'] = 2,  ['D'] = 3,  ['E'] = 4,  ['F'] = 5,
//    ['G'] = 6,  ['H'] = 7,  ['I'] = 8,  ['K'] = 9,  ['L'] = 10,
//    ['M'] = 11, ['N'] = 12, ['P'] = 13, ['Q'] = 14, ['R'] = 15,
//    ['S'] = 16, ['T'] = 17, ['V'] = 18, ['W'] = 19, ['Y'] = 20,
//    ['-'] = 0
//};

float get_ctilde(const vector<string>& sequences, int i, int j, int Q) {
    vector<int> onept_count_i(Q, 0);
    vector<int> onept_count_j(Q, 0);
    vector<vector<int>> twopt_counts(Q, vector<int>(Q, 0));
    
    const int n_seqs = sequences.size();
    
    // Count occurrences
    for (const auto& seq : sequences) {
        int i_aa_idx = AA_LOOKUP[seq[i]];
        int j_aa_idx = AA_LOOKUP[seq[j]];
        
        onept_count_i[i_aa_idx]++;
        onept_count_j[j_aa_idx]++;
        twopt_counts[i_aa_idx][j_aa_idx]++;
    }
    
    // Calculate correlation
    float ctilde = 0.0f;
    const float n_seqs_f = static_cast<float>(n_seqs);
    
    for (int k = 0; k < Q; ++k) {
        for (int l = 0; l < Q; ++l) {
            float freq_pair = twopt_counts[k][l] / n_seqs_f;
            float freq_prod = (onept_count_i[k] * onept_count_j[l]) / (n_seqs_f * n_seqs_f);
            float diff = freq_pair - freq_prod;
            ctilde += diff * diff;
        }
    }
    
    return sqrt(ctilde);
}

// Pre-read the alignment file into memory
vector<string> read_alignment(const string& filename) {
    vector<string> sequences;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return sequences;
    }
    
    string line;
    while (getline(file, line)) {
        sequences.push_back(line);
    }
    file.close();
    return sequences;
}

int main() {
    const int N = 262;  // Matrix dimension
    const int Q = 21;   // Number of amino acids + gap
    
    // Pre-read the alignment
    vector<string> sequences = read_alignment("iter_aln_short.txt");
    if (sequences.empty()) return 1;
    
    // Allocate result matrix
    vector<vector<float>> ctilde_arr(N, vector<float>(N));
    
    // Parallel processing of matrix elements
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            ctilde_arr[i][j] = get_ctilde(sequences, i, j, Q);
            
            // Print progress only at multiples of 100 and in the master thread
            if (omp_get_thread_num() == 0 && i % 100 == 0 && j % 100 == 0) {
                cout << "Progress: i=" << i << ", j=" << j << endl;
            }
        }
    }
    
    // Print final progress
    cout << "Completed: i=" << N-1 << ", j=" << N-1 << endl;
    
    // Write results to file
    ofstream ctilde_file("ctilde_iter.csv");
    if (!ctilde_file.is_open()) {
        cerr << "Error opening output file!" << endl;
        return 1;
    }
    
    ctilde_file << fixed << setprecision(6);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            ctilde_file << ctilde_arr[i][j] << ",";
        }
        ctilde_file << ctilde_arr[i][N-1] << "\n";
    }
    
    return 0;
}
