#include "cords.hpp"
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <random>
#include <chrono>
#include <functional>
#include "units_pair.hpp"

namespace algos {

    std::vector<std::vector<std::string>> Cords::GetTable(const std::string &filename) {
        std::ifstream input;
        input.open(filename.c_str());
        if (input.is_open()) {
            std::string line;
            std::vector<std::vector<std::string>> res;
            while (getline(input, line)) {
                std::vector<std::string> tmp;
                std::string word;
                for (char c: line) {
                    if (c != '|') {
                        word.push_back(c);
                    } else {
                        tmp.push_back(word);
                        word.clear();
                    }
                }
                tmp.push_back(word);
                res.push_back(tmp);
            }
            input.close();
            return res;
        }
        std::cout << "couldn't open" << std::endl;
        return {};

    }

    Cords::Cords() = default;

    //TODO: скорее всего можно вынести в парсинг
    unsigned long long Cords::ComputeCardinality(size_t ind, const Table &data) {
        std::unordered_map<std::string, bool> map; //TODO: понять какой тип будет когда буду мержить
        unsigned long long res = 0;
        for (size_t i = 1; i < data.size(); i++) {
            if (map.find(data[i][ind]) == map.end()) {
                map.emplace(data[i][ind], true);
                res++;
            }
        }
        //std::cout << "cardinality of " << data[0][ind] << " is " << res << std::endl;
        return res;
    }

    //TODO: computes correctly only for 0<p<0.4
    //TODO: попробовать убрать unsigned
// formulae (2) from "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies."
    unsigned long long Cords::ComputeSampleSize(unsigned long long d1, unsigned long long d2) const {
        long double v = (d1 - 1) * (d2 - 1);
        long double d = std::min(d1, d2);
        long double log = std::log(p * std::sqrt(2 * pi));
        long double numerator = std::pow(-16 * v * log, 0.5) - 8 * log;
        long double denumerator = delta * (d - 1);
        long double v2 = std::pow(v, 0.071);
        return (numerator / denumerator) * (v2 / 1.69);
        //return ((std::pow(-16 * v * log, 0.5) - 8 * log) / (1.69 * delta * (d - 1) * std::pow(v, -0.071)));
    }

    Sample Cords::Sampling(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table &data) {
        Sample res;
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        auto dice_rand = std::bind(std::uniform_int_distribution<size_t>(1, rows - 1),
                                   std::mt19937(seed));
        std::unordered_map<std::string, bool> map_lhs;
        std::unordered_map<std::string, bool> map_rhs;
        size_t card_lhs = 0;
        size_t card_rhs = 0;
        for (int i = 0; i < sample_size; i++) {
            size_t ind = dice_rand();
            res.self.emplace_back(l, r, ind);
            if (map_lhs.find(data[i][ind]) == map.end()) {
                map.emplace(data[i][ind], true);
                res++;
            }
        }
        return res;
    }

    void PrintSample(const std::vector<UnitsPair> &sample, const Cords::Table &tbl) {
        for (auto i: sample) {
            std::cout << tbl[i.row_ind][i.col_l] << " " << tbl[i.row_ind][i.col_r] << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    //TODO: Понять нужно ли сортировать столбцы по количеству различных значений, скорее всего надо, но лучше использовать сортировку которая за один проход поймёт что ничего делать не надо
    //TODO: Думаю надо модифицировать квик сорт либо которая ведёт элемент назад пока не встанет куда надо
    void Cords::ExecuteInternal(const std::string &filename) {
        if (eps1 >= delta) {
            throw std::logic_error("delta must be more than eps1");
        }
        auto data = GetTable(filename);

        for (const auto &i: data[0]) {
            std::cout << i << " ";
        }

        std::cout << std::endl;
        size_t rows = data.size();
        size_t columns = data[0].size();
        std::cout << "columns: " << columns << " rows: " << rows << std::endl;
        std::vector<bool> is_soft(columns, false);
        std::vector<unsigned long long> cardinality(columns, 0);
        // every pair of columns in which one of the columns is almost a key //TODO: not sure if the description is correct
        for (size_t col_i = 0; col_i < columns - 1; col_i++) {
            if (is_soft[col_i]) {
                continue;
            }
            if (cardinality[col_i] == 0) {
                cardinality[col_i] = ComputeCardinality(col_i, data);
            }
            if (cardinality[col_i] >= (1 - eps1) * rows) {
                soft_keys.push_back(col_i);
                is_soft[col_i] = true;
                continue;
            }
            //TODO: check whether column col_i is trivial

            for (size_t col_k = col_i + 1; col_k < columns; col_k++) {
                if (is_soft[col_k]) {
                    continue;
                }
                if (cardinality[col_k] == 0) {
                    cardinality[col_k] = ComputeCardinality(col_k, data);
                }
                if (cardinality[col_k] >= (1 - eps1) * rows) {
                    soft_keys.push_back(col_k);
                    is_soft[col_k] = true;
                    //TODO: check whether column col_k is trivial
                    continue;
                }

                unsigned long long sample_size = ComputeSampleSize(cardinality[col_i], cardinality[col_k]);
                size_t size = rows > sample_size ? sample_size : rows;
                auto sample = Sampling(size, rows, col_i, col_k, data);
                //TODO: скорее всего надо будет передавать не rows а минимум из длины
                //PrintSample(sample, data);
            }
        }
//        std::cout << std::endl;
//        for (int i = 0; i < columns; i++) {
//            std::cout << "cardinality for " << data[0][i] << " is: " << cardinality[i] << std::endl;
//        }
        std::cout << "soft keys are: ";
        for (auto i: soft_keys) {
            std::cout << data[0][i] << " ";
        }
    }
}