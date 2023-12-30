#pragma once

#include <vector>
#include <random>
#include <chrono>
#include <functional>
#include <iostream>
#include <unordered_set>

namespace algos {
    using Table = std::vector<std::vector<std::string>>;

    struct Sample {
        std::vector<size_t> row_indices;
        size_t lhs_ind;
        size_t rhs_ind;
        size_t lhs_cardinality;
        size_t rhs_cardinality;
        size_t concat_cardinality;

        Sample(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table &data, bool has_header) {
            lhs_ind = l;
            rhs_ind = r;

            auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            auto dice_rand = std::bind(std::uniform_int_distribution<size_t>(has_header ? 1 : 0, rows - 1),
                                       std::mt19937(seed));

            std::unordered_set<std::string> map_lhs; //using unordered_set because we need average const time complexity for find and emplace
            std::unordered_set<std::string> map_rhs;
            std::unordered_set<std::string> map_cardinality;


            for (int i = 0; i < sample_size; i++) {
                size_t row = dice_rand(); // номер строки!!!!
                row_indices.push_back(row);
                if (map_lhs.find(data[row][l]) ==
                    map_lhs.end()) { //computes cardinality of a left-hand column in a sample
                    map_lhs.emplace(data[row][l]);
                }
                if (map_rhs.find(data[row][r]) ==
                    map_rhs.end()) { //computes cardinality of a right-hand column in a sample
                    map_rhs.emplace(data[row][r]);
                }

                //TODO: might need fixes when integrated to desbordante due to custom type system
                std::string tmp = data[row][l] + data[row][r];
                if (map_cardinality.find(tmp) ==
                    map_cardinality.end()) {// computes cardinality of columns concatenation
                    map_cardinality.emplace(tmp);
                }
            }
            lhs_cardinality = map_lhs.size();
            rhs_cardinality = map_rhs.size();
            concat_cardinality = map_cardinality.size();
        }
    };
}