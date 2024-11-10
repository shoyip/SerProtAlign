#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <functional>

// Custom hash function for unordered_set<int>
struct set_hash {
    template <typename T>
    size_t operator()(const std::unordered_set<T>& set) const {
        size_t seed = 0;
        for (const auto& elem : set) {
            seed ^= std::hash<T>{}(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

int main() {
    std::vector<std::vector<char>> transposed_alignment;

    std::ifstream file("aln_short.txt");

    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(file, line)) {
        transposed_alignment.resize(line.size());
        std::vector<char> sequence(line.begin(), line.end());
        for (size_t col=0; col<line.size(); ++col) {
            transposed_alignment[col].push_back(sequence[col]);
        }
    }

    file.close();

//    std::cout << "The transposed 2D alignment is:" << std::endl;
//    for (const auto& column : transposed_alignment) {
//        for (const auto& residue : column) {
//            std::cout << residue;
//        }
//        std::cout << std::endl;
//    }
//
//
//    for (size_t j=0; j<transposed_alignment.size(); ++j) {
//        for (size_t i=0; i<transposed_alignment[j].size(); ++i) {
//            std::cout << "i" << i << "j" << j << "res" << transposed_alignment[j][i] << std::endl;
//        }
//    }

    std::vector<std::map<char, int>> char_counts_alignment;
    for (size_t j=0; j<transposed_alignment.size(); ++j) {
        std::map<char, int> char_counts;
        for (size_t i=0; i<transposed_alignment[j].size(); ++i) {
            char_counts[transposed_alignment[j][i]] += 1;
        }
        char_counts_alignment.push_back(char_counts);
//        for (const auto& [key, value]: char_counts) {
//            std::cout  << "Character " << key << " has frequency " << value * 1. / 5. << std::endl;
//        }
//        std::cout << "End of column." << std::endl;
    }

    std::vector<std::vector<std::unordered_map<std::unordered_set<char>,int,set_hash>>> twopoint_counts;
    for (size_t i=0; i<transposed_alignment.size(); ++i) {
//    for (const auto& column_i : transposed_alignment) {
        std::vector<std::unordered_map<std::unordered_set<char>,int,set_hash>> twopoint_counts_row;
        for (size_t j=0; j<transposed_alignment.size(); ++j) {
//        for (const auto& column_j : transposed_alignment) {        
            std::unordered_map<std::unordered_set<char>,int,set_hash> twopoint_count;
            for (size_t k=0; k<transposed_alignment[i].size(); ++k) {
                //since it's a set, it will show just one letter if these are identical i.e. II->I
                std::unordered_set<char> char_tuple = {transposed_alignment[i][k], transposed_alignment[j][k]};
                twopoint_count[char_tuple] += 1;
            }
            twopoint_counts_row.push_back(twopoint_count);
        }
        twopoint_counts.push_back(twopoint_counts_row);
    }

    for (const auto& [key, value] : twopoint_counts[1][2]) {
        for (const auto& residue : key) {
//            int aa_idx;
//            switch (residue) {
//                case '-':
//                    aa_idx = 0;
//                    break;
//                case 'A':
//                    aa_idx = 1;
//                    break;
//                case 'C':
//                    aa_idx = 2;
//                    break;
//                case 'D':
//                    aa_idx = 3;
//                    break;
//                case 'E':
//                    aa_idx = 4;
//                    break;
//                case 'F':
//                    aa_idx = 5;
//                    break;
//                case 'G':
//                    aa_idx = 6;
//                    break;
//                case 'H':
//                    aa_idx = 7;
//                    break;
//                case 'I':
//                    aa_idx = 8;
//                    break;
//                case 'K':
//                    aa_idx = 9;
//                    break;
//                case 'L':
//                    aa_idx = 10;
//                    break;
//                case 'M':
//                    aa_idx = 11;
//                    break;
//                case 'N':
//                    aa_idx = 12;
//                    break;
//                case 'P':
//                    aa_idx = 13;
//                    break;
//                case 'Q':
//                    aa_idx = 14;
//                    break;
//                case 'R':
//                    aa_idx = 15;
//                    break;
//                case 'S':
//                    aa_idx = 16;
//                    break;
//                case 'T':
//                    aa_idx = 17;
//                    break;
//                case 'V':
//                    aa_idx = 18;
//                    break;
//                case 'W':
//                    aa_idx = 19;
//                    break;
//                case 'Y':
//                    aa_idx = 20;
//                    break;
//            }
//            std::cout << aa_idx << std::endl;
        }
        std::cout << value << std::endl;
    }

    return 0;
}
