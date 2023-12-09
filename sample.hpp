#pragma once

#include <vector>
#include <random>
#include <chrono>
#include <functional>
#include <iostream>
#include "alliacies.hpp"
#include "unordered_set"

namespace algos {
    struct Sample {
        std::vector<size_t> row_indices;
        size_t lhs_ind;
        size_t rhs_ind;
        size_t lhs_cardinality;
        size_t rhs_cardinality;
        size_t concat_cardinality;

        Sample(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table &data) {
            lhs_ind = l;
            rhs_ind = r;

            auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            auto dice_rand = std::bind(std::uniform_int_distribution<size_t>(1, rows - 1),
                                       std::mt19937(seed));
            std::unordered_set<std::string> map_lhs;
            std::unordered_set<std::string> map_rhs;
            std::unordered_set<std::string> map_cardinality;
            size_t card_lhs = 0;
            size_t card_rhs = 0;
            size_t card_concat = 0;
            for (int i = 0; i < sample_size; i++) {
                size_t ind = dice_rand(); //- номер строки !!!!
                row_indices.push_back(ind);
                if (map_lhs.find(data[ind][l]) == map_lhs.end()) {
                    map_lhs.emplace(data[ind][l]);
                    card_lhs++;
                }
                if (map_rhs.find(data[ind][r]) == map_rhs.end()) {
                    map_rhs.emplace(data[ind][r]);
                    card_rhs++;
                }
                std::string tmp = data[ind][l] + data[ind][r];
                if (map_cardinality.find(tmp) == map_cardinality.end()) {
                    map_cardinality.emplace(tmp);
                    card_concat++;
                }
            }

            lhs_cardinality = card_lhs;
            rhs_cardinality = card_rhs;
            concat_cardinality = card_concat;
        }
    };
}