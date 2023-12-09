#include "cords.hpp"
#include <iostream>
#include <unordered_map>
#include <cmath>
#include "sample.hpp"
#include <algorithm>
#include <functional>
#include <boost/math/distributions/chi_squared.hpp>

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

    //TODO: computes correctly only for 0<p<0.4
    //TODO: попробовать убрать unsigned
// formulae (2) from "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies."
    long long Cords::ComputeSampleSize(unsigned long long d1, unsigned long long d2) const {
        long double v = (d1 - 1) * (d2 - 1);
        long double D = std::min(d1, d2);
        long double log = std::log(p * std::sqrt(2 * pi));
        long double numerator = std::pow(-16 * v * log, 0.5) - 8 * log;
        long double denumerator = delta * (D - 1);
        long double v2 = std::pow(v, 0.071);
        return static_cast<long long>((numerator / denumerator) * (v2 / 1.69));
        //return ((std::pow(-16 * v * log, 0.5) - 8 * log) / (1.69 * delta * (d - 1) * std::pow(v, -0.071)));
    }

    void PrintSample(const Sample &sample, const Table &tbl) {
        for (auto i: sample.row_indices) {
            std::cout << tbl[i][sample.lhs_ind] << " " << tbl[i][sample.rhs_ind] << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    bool Cords::DetectSfd(const Sample &smp) {
        if (smp.concat_cardinality <= eps2 * smp.row_indices.size()) {
            if (smp.lhs_cardinality >= (1 - eps3) * smp.concat_cardinality) {
                sfd.emplace_back(smp.lhs_ind, smp.rhs_ind);
                return true;
            } else if (smp.rhs_cardinality >= (1 - eps3) * smp.concat_cardinality) {
                sfd.emplace_back(smp.rhs_ind, smp.lhs_ind);
                return true;
            }
        }
        return false;
    }

    void Cords::GetFrequentValuesStatistics(const Table &data) {
        for (size_t col_ind = 0; col_ind < data[0].size(); col_ind++) {
            std::unordered_map<std::string, int> counter;
            //TODO: CHANGE
            size_t st = has_header ? 1 : 0;
            for (size_t row_ind = st; row_ind < data.size(); row_ind++) {
                if (counter.find(data[row_ind][col_ind]) == counter.end()) {
                    counter[data[row_ind][col_ind]] = 1;
                    cardinality[col_ind]++;
                } else {
                    counter[data[row_ind][col_ind]]++;
                }
            }
            //TODO: CHANGE!!!!!!!!!!!!!!!!
            std::vector<std::pair<std::string, int>> stats;
            for (auto &it: counter) {
                stats.emplace_back(it.first, it.second);
            }

            auto cmp = [](const std::pair<std::string, int> &left, const std::pair<std::string, int> &right) {
                return left.second > right.second;
            };
            std::sort(stats.begin(), stats.end(), cmp);

            for (int i = 0; i < std::min(num_freqvals, stats.size()); i++) {
                frequent_values_statistics[col_ind][stats[i].first] = {stats[i].second, i};
                freq_sums[col_ind] += stats[i].second;
            }
        }
    }

    void
    Cords::Filter(const std::unordered_map<std::string, std::pair<size_t, size_t>> &freq_val_stats,
                  const algos::Table &data, size_t col_ind,
                  Sample &smp) {
        for (auto row = smp.row_indices.begin(); row != smp.row_indices.end(); row++) {
            if (freq_val_stats.find(data[*row][col_ind]) == freq_val_stats.end()) {
                smp.row_indices.erase(row);
                std::cout << "erased: " << data[*row][col_ind] << std::endl;
            }
        }

    }

    size_t Cords::Category(size_t col_ind, const std::string &val, size_t domain, bool skew, const Table &data) const {
        if (skew) {
            return frequent_values_statistics[col_ind].at(val).second;
        }
        //TODO: CHANGE HASH
        return std::hash<std::string>{}(val) % domain; //+1;
    }

    size_t IsZero(size_t val) {
        return val == 0;
    }

    //TODO: Понять нужно ли сортировать столбцы по количеству различных значений, скорее всего надо, но лучше использовать сортировку которая за один проход поймёт что ничего делать не надо
    //TODO: Думаю надо модифицировать квик сорт либо которая ведёт элемент назад пока не встанет куда надо
    void Cords::ExecuteInternal(const std::string &filename, bool has_head) {
        this->has_header = has_head;
        int cnt = 0;
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
        cardinality.resize(columns, 0);
        frequent_values_statistics.resize(columns);
        freq_sums.resize(columns, 0);
        is_skewed.resize(columns, false);
        domains.resize(columns, 0);
        GetFrequentValuesStatistics(data);
        for (size_t col_i = 0; col_i < columns - 1; col_i++) {
            if (is_soft[col_i]) {
                continue;
            }
            if (cardinality[col_i] >= (1 - eps1) * rows) {
                soft_keys.push_back(col_i);
                is_soft[col_i] = true;
                continue;
            }
            if (cardinality[col_i] == 1) {
                trivial_columns.push_back(col_i);
                continue;
            }

            for (size_t col_k = col_i + 1; col_k < columns; col_k++) {
                if (is_soft[col_k]) {
                    continue;
                }
                if (cardinality[col_k] >= (1 - eps1) * rows) {
                    soft_keys.push_back(col_k);
                    is_soft[col_k] = true;
                    continue;
                }
                if (cardinality[col_k] == 1) {
                    trivial_columns.push_back(col_k);
                    continue;
                }
                cnt++;
                unsigned long long sample_size = ComputeSampleSize(cardinality[col_i], cardinality[col_k]);

                size_t size = rows > sample_size ? sample_size : rows;
                Sample smp(size, rows, col_i, col_k, data);
                //TODO: скорее всего надо будет передавать не rows а минимум из длины
                //PrintSample(smp, data);

                if (DetectSfd(smp)) {
                    continue;
                }

                //skew handling
                for (size_t col_ind: {col_i, col_k}) {
                    if (freq_sums[col_ind] >= (1 - eps4) * data.size()) {
                        is_skewed[col_ind] = true;
                        domains[col_ind] = frequent_values_statistics[col_ind].size();
                        Filter(frequent_values_statistics[col_ind], data, col_ind, smp);
                    } else {
                        domains[col_ind] = std::min(cardinality[col_ind], d_max);
                    }
                }
                //std::cout << "skew handled" << std::endl;

                std::vector<std::vector<size_t>> n_i_j(domains[col_i], std::vector<size_t>(domains[col_k], 0));
                std::vector<size_t> n_i(domains[col_i], 0);
                std::vector<size_t> n_j(domains[col_k], 0);

                //std::cout << "d1: " << domains[col_i] << " d2: " << domains[col_k] << std::endl;
                //std::cout << col_i << ":" << is_skewed[col_i] << " " << col_k << ":" << is_skewed[col_k] << std::endl;
                for (size_t row_ind: smp.row_indices) {
                    size_t i = Category(col_i, data[row_ind][col_i], domains[col_i], is_skewed[col_i], data);
                    size_t j = Category(col_k, data[row_ind][col_k], domains[col_k], is_skewed[col_k], data);
                    //std::cout << "i: " << i << " j: " << j << std::endl;
                    n_i_j[i][j]++;
                    n_i[i]++;
                    n_j[j]++;
                }

//                std::cout << "n_i_j:" << std::endl;
//                for (const auto &it: n_i_j) {
//                    std::cout << "{";
//                    for (auto it2: it) {
//                        std::cout << it2 << " ";
//                    }
//                    std::cout << "}" << std::endl;
//                }
//                std::cout << "n_i:" << std::endl;
//                std::cout << "{";
//                for (auto it: n_i) {
//                    std::cout << it << " ";
//                }
//                std::cout << "}" << std::endl;


                size_t sum = 0;
                for (int i = 0; i < domains[col_i]; i++) {
                    for (int j = 0; j < domains[col_k]; j++) {
                        sum += IsZero(n_i_j[i][j]);
                    }
                }
                std::cout << "zeros sum:" << sum << std::endl;
                if (sum > eps5 * domains[col_i] * domains[col_k]) {
                    correlations.emplace_back(col_i, col_k);
                    //std::cout << "found correlation on eps: " << col_i << " <->" << col_k << std::endl;
                    continue;
                }

                long double chi_squared = 0;
                for (int i = 0; i < domains[col_i]; i++) {
                    for (int j = 0; j < domains[col_k]; j++) {
//                        std::cout << "n_i_j:" << n_i_j[i][j] << std::endl;
//                        std::cout << "n_i:" << n_i[i] << std::endl;
//                        std::cout << "n_j:" << n_j[j] << std::endl;
//                        std::cout << "res: " << ((n_i_j[i][j] - n_i[i] * n_j[j]) * (n_i_j[i][j] - n_i[i] * n_j[j])) /
//                                                (n_i[i] * n_j[j]) << std::endl;
                        std::cout << "\n\n" << std::endl;
                        chi_squared += // formulae (1) from "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies."
                                ((n_i_j[i][j] - n_i[i] * n_j[j]) * (n_i_j[i][j] - n_i[i] * n_j[j])) / (n_i[i] * n_j[j]);
                    }
                }

                size_t v = (domains[col_i] - 1) * (domains[col_k] - 1);
                std::cout << "degrees of freedom:" << v << std::endl;
                boost::math::chi_squared dist(v);
                long double t = quantile(dist, 1 - p);
                //std::cout << "chi squared: " << chi_squared << " quantile: " << t << std::endl;
                if (chi_squared > t) {
                    correlations.emplace_back(col_i, col_k);
                    continue;
                }
            }
        }


        std::cout << "soft keys are: ";
        for (auto i: soft_keys) {
            std::cout << data[0][i] << " ";
        }
        std::cout << std::endl;
        std::cout << "SFD ARE" << std::endl;
        for (auto i: sfd) {
            std::cout << data[0][i.first] << "=>" << data[0][i.second] << std::endl;
        }
        std::cout << "overall sfd:" << sfd.size() << std::endl;
        std::cout << "overall pairs:" << cnt << std::endl;
        std::cout << "correlations: " << std::endl;
        for (auto it: correlations) {
            std::cout << "{" << it.first << "," << it.second << "}" << " ";
        }
        std::cout << std::endl;

    }
}