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
#include <random>
#include <algorithm>
#include <unordered_set>

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

vector<int> generate_random_indices(int list_length, int max_int, unsigned int seed = std::random_device{}()) {
    // Create a random number generator
    mt19937 gen(seed);  // Mersenne Twister random number engine

    // Create a vector of all possible indices
    vector<int> all_indices(max_int);
    iota(all_indices.begin(), all_indices.end(), 1);  // Fill with 1, 2, 3, ..., N

    // Shuffle the indices
    shuffle(all_indices.begin(), all_indices.end(), gen);

    // Take the first 'listLength' elements
    return vector<int>(all_indices.begin(), all_indices.begin() + list_length);
}

template<typename T>
vector<T> select_elements(const vector<T>& source_vector, const vector<int>& indices) {
    vector<T> selected_elements;

    for (int index : indices) {
        // Adjust index to 0-based indexing and handle out-of-bounds
        if (index > 0 && index <= source_vector.size()) {
            selected_elements.push_back(source_vector[index - 1]);
        }
    }

    return selected_elements;
}

float get_ctilde(const vector<string>& sequences, const vector<float>& weights, int i, int j, int Q) {
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
    
    return weights[i] * weights[j] * sqrt(ctilde);
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
        if (!line.empty() && line[0] != '>') {
            sequences.push_back(line);
        }
    }
    file.close();
    return sequences;
}

void shuffle_alignment(vector<string>& alignment, unsigned int seed = std::random_device{}()) {
    if (alignment.empty() || alignment[0].empty()) {
        throw invalid_argument("Sequences must not be empty.");
    }


    size_t num_columns = alignment[0].size();
    for (const auto& seq : alignment) {
        if (seq.size() != num_columns) {
            throw invalid_argument("All sequences must be the same length.");
        }
    }

    vector<string> columns(num_columns);
    for (size_t col = 0; col < num_columns; ++col) {
        for (const auto& seq : alignment) {
            columns[col] += seq[col];
        }
    }

    mt19937 rng(seed);
    for (auto& column : columns) {
        shuffle(column.begin(), column.end(), rng);
    }

    for (size_t row = 0; row < alignment.size(); ++row) {
        for (size_t col = 0; col < num_columns; ++col) {
            alignment[row][col] = columns[col][row];
        }
    }
}

vector<float> read_weights(const string& filename) {
    vector<float> weights;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return weights;
    }
    float weight;
    while (file >> weight) {
        weights.push_back(weight);
    }
    file.close();
    return weights;
}

int main(int argc, char* argv[]) {
    const int N = atoi(argv[1]);  // Matrix dimension
    const int Q = 21;   // Number of amino acids + gap

    // Pre-read the alignment
    vector<string> full_sequences = read_alignment(argv[2]);
    if (full_sequences.empty()) return 1;

    vector<float> full_weights = read_weights(argv[3]);
    if (full_weights.empty()) {
        full_weights.resize(full_sequences.size(), 1);
    }
    
    vector<int> random_indices = generate_random_indices(atoi(argv[4]), full_sequences.size(), 42);
    vector<string> sequences = select_elements(full_sequences, random_indices);
    vector<float> weights = select_elements(full_weights, random_indices);
    
    // Allocate result matrix
    vector<vector<float>> ctilde_arr(N, vector<float>(N));
    
    // Parallel processing of matrix elements
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            ctilde_arr[i][j] = get_ctilde(sequences, weights, i, j, Q);
            
            // Print progress only at multiples of 100 and in the master thread
            if (omp_get_thread_num() == 0 && i % 100 == 0 && j % 100 == 0) {
                cout << "Progress: i=" << i << ", j=" << j << endl;
            }
        }
    }
    
    // Print final progress
    cout << "Completed: i=" << N-1 << ", j=" << N-1 << endl;
    
    // Write results to file
    ofstream ctilde_file(argv[5]);
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

    cout << "Starting computation of statistics for shuffled alignment." << endl;

    shuffle_alignment(sequences, 42);

    // Ctilde for shuffled alignment
    vector<vector<float>> ctilde_shuffled_arr(N, vector<float>(N));
    
    // Parallel processing of matrix elements
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            ctilde_shuffled_arr[i][j] = get_ctilde(sequences, weights, i, j, Q);
            
            // Print progress only at multiples of 100 and in the master thread
            if (omp_get_thread_num() == 0 && i % 100 == 0 && j % 100 == 0) {
                cout << "Progress: i=" << i << ", j=" << j << endl;
            }
        }
    }
    
    // Print final progress
    cout << "Completed: i=" << N-1 << ", j=" << N-1 << endl;
    
    // Write results to file
    ofstream ctilde_shuffled_file(argv[6]);
    if (!ctilde_shuffled_file.is_open()) {
        cerr << "Error opening output file!" << endl;
        return 1;
    }
    
    ctilde_shuffled_file << fixed << setprecision(6);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            ctilde_shuffled_file << ctilde_shuffled_arr[i][j] << ",";
        }
        ctilde_shuffled_file << ctilde_shuffled_arr[i][N-1] << "\n";
    }
    
    return 0;
}
