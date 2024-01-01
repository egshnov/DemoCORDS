#include <iostream>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <boost/math/distributions/chi_squared.hpp>

#include "sample.hpp"
#include "cords.hpp"

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

    // formulae (2) from "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies."
    //FIXME: computes correctly only for 0<p<0.4
    long long Cords::ComputeSampleSize(long long d1, long long d2) const {
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
        std::cout << "sample: " << std::endl;
        for (auto i: sample.row_indices) {
            std::cout << tbl[i][sample.lhs_ind] << " " << tbl[i][sample.rhs_ind] << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    void PrintContingencyTable(const std::vector<std::vector<long double>> &table) {
        for (const auto &it: table) {
            for (auto it2: it) {
                std::cout << it2 << " ";
            }
            std::cout << std::endl;
        }
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

    //calculates cardinality of each column in a table, creates a vector of unordered_maps comprised of num_freqvals most frequent values
    // where key = value val = occurrence frequency
    void Cords::GetFrequentValuesStatistics(const Table &data) {
        for (size_t col_ind = 0; col_ind < data[0].size(); col_ind++) {
            std::unordered_map<std::string, int> counter;
            //TODO: CHANGE has header?
            size_t st = has_header ? 1 : 0;
            for (size_t row_ind = st; row_ind < data.size(); row_ind++) {
                auto it = counter.find(data[row_ind][col_ind]);
                if (it == counter.end()) {
                    counter[data[row_ind][col_ind]] = 1;
                    cardinality[col_ind]++;
                } else {
                    it->second++;
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
                frequent_values_statistics[col_ind][stats[i].first] = {stats[i].second,
                                                                       i}; // key = value map[key] = {value frequency, sequence number in the column}
                freq_sums[col_ind] += stats[i].second;
            }
        }
    }

    void Cords::Filter(const std::unordered_map<std::string, std::pair<size_t, size_t>> &freq_val_stats,
                       const algos::Table &data, size_t col_ind,
                       Sample &smp) {
        auto row = smp.row_indices.begin();
        while (row != smp.row_indices.end()) {
            if (freq_val_stats.find(data[*row][col_ind]) == freq_val_stats.end()) {
                std::cout << "erased: " << *row << " " << data[*row][col_ind] << std::endl;
                row = smp.row_indices.erase(row);
            } else {
                row++;
            }
        }
    }

    //FIXME: Redundant changes in domains[col_i] and is_skewed[col_i] in case if we already know if column is skewed or not
    void Cords::SkewHandling(size_t col_i, size_t col_k, const Table &data, Sample &smp) {
        for (size_t col_ind: {col_i, col_k}) {
            if (freq_sums[col_ind] >= (1 - eps4) * data.size()) {
                is_skewed[col_ind] = true;
                domains[col_ind] = frequent_values_statistics[col_ind].size();
                Filter(frequent_values_statistics[col_ind], data, col_ind, smp);
            } else {
                domains[col_ind] = std::min(cardinality[col_ind], d_max);
            }
        }
    }

    size_t Cords::Category(size_t col_ind, const std::string &val, size_t domain, bool skew, const Table &data) const {
        if (skew) {
            return frequent_values_statistics[col_ind].at(val).second;
        }
        //TODO: change hash when adding to desbordante
        return std::hash<std::string>{}(val) % domain; //+1;
    }

    void Cords::ComputeContingencyTable(Sample &smp, const Table &data, size_t col_i, size_t col_k,
                                        std::vector<std::vector<long double>> &n_i_j, std::vector<long double> &n_i,
                                        std::vector<long double> &n_j) {
        for (size_t row_ind: smp.row_indices) {
            size_t i = Category(col_i, data[row_ind][col_i], domains[col_i], is_skewed[col_i], data);
            size_t j = Category(col_k, data[row_ind][col_k], domains[col_k], is_skewed[col_k], data);
            n_i_j[i][j]++;
            n_i[i]++;
            n_j[j]++;
        }
    }
    //TODO: change val to long double?
    size_t IsZero(size_t val) {
        return val == 0;
    }

    bool Cords::TooMuchStructuralZeros(const std::vector<std::vector<long double>> &n_i_j, size_t col_i, size_t col_k) {
        long double zeros_sum = 0;
        for (int i = 0; i < domains[col_i]; i++) {
            for (int j = 0; j < domains[col_k]; j++) {
                zeros_sum += IsZero(n_i_j[i][j]);
            }
        }
        return zeros_sum > eps5 * domains[col_i] * domains[col_k];
    }

    long double
    Cords::ComputeChiSquared(const std::vector<std::vector<long double>> &n_i_j, const std::vector<long double> &n_i,
                             const std::vector<long double> &n_j, size_t col_i, size_t col_k, long double sample_size) {
        long double chi_squared = 0;
        for (size_t i = 0; i < domains[col_i]; i++) {
            for (size_t j = 0; j < domains[col_k]; j++) {
                if (n_i[i] * n_j[j] == 0) throw std::logic_error("Structural zeroes are not handled");
                long double actual = n_i_j[i][j];
                long double expected = n_i[i] * n_j[j] / sample_size;
                chi_squared += (actual - expected) * (actual - expected) / (expected);
                // similar to formulae (1) from "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies."
                // seems like formulae which is given in the paper is incorrect
            }
        }
        return chi_squared;
    }

    void Cords::Init(size_t columns) {
        cardinality.resize(columns, 0);
        frequent_values_statistics.resize(columns);
        freq_sums.resize(columns, 0);
        is_skewed.resize(columns, false);
        domains.resize(columns, 0);
    }

    std::pair<size_t, size_t> Cords::SwitchIndicesIfNecessary(size_t ind1, size_t ind2) {
        if (cardinality[ind2] > cardinality[ind1]) {
            return {ind2, ind1};
        }
        return {ind1, ind2};
    }

    void Cords::ExecuteInternal(const std::string &filename, bool has_head) {
        has_header = has_head;
        int cnt = 0;
        if (eps1 >= delta) {
            throw std::logic_error("delta must be more than eps1");
        }
        auto data = GetTable(filename);
        if (has_header) {
            for (const auto &i: data[0]) {
                std::cout << i << " ";
            }
        } else {
            std::cout << "NO HEADER" << std::endl;
        }

        std::cout << std::endl;

        size_t rows = data.size();
        size_t columns = data[0].size();
        std::vector<bool> is_soft(columns, false);
        std::cout << "columns: " << columns << " rows: " << rows << std::endl;

        Init(columns);

        GetFrequentValuesStatistics(data);

        for (auto it: cardinality) {
            std::cout << it << " ";
        }
        std::cout << std::endl;
        for (size_t ind1 = 0; ind1 < columns - 1; ind1++) {

            if (is_soft[ind1]) {
                continue;
            }
            if (cardinality[ind1] >= (1 - eps1) * rows) {
                soft_keys.push_back(ind1);
                is_soft[ind1] = true;
                continue;
            }
            if (cardinality[ind1] == 1) {
                trivial_columns.push_back(ind1);
                continue;
            }

            for (size_t ind2 = ind1 + 1; ind2 < columns; ind2++) {

                if (is_soft[ind2]) {
                    continue;
                }
                if (cardinality[ind2] >= (1 - eps1) * rows) {
                    soft_keys.push_back(ind2);
                    is_soft[ind2] = true;
                    continue;
                }
                if (cardinality[ind2] == 1) {
                    trivial_columns.push_back(ind2);
                    continue;
                }

                auto [col_i, col_k] = SwitchIndicesIfNecessary(ind1, ind2);
                cnt++;
                unsigned long long possible_sample_size = ComputeSampleSize(cardinality[col_i], cardinality[col_k]);


                std::cout << "possible sample size: " << possible_sample_size << std::endl;
                //size_t size = rows > possible_sample_size ? possible_sample_size : rows;
                //  Sample smp(size, rows, col_i, col_k, data);

                Sample smp(possible_sample_size, rows, col_i, col_k, data, has_header);
                //TODO: скорее всего надо будет передавать не rows а минимум из длинн upd: похоже что нет, подумать ещё
                //PrintSample(smp, data);

                if (DetectSfd(smp)) {
                    continue;
                }
                SkewHandling(col_i, col_k, data, smp);
                //std::cout << "filtered ";
                //PrintSample(smp, data);
                std::vector<std::vector<long double>> n_i_j(domains[col_i],
                                                            std::vector<long double>(domains[col_k], 0));
                std::vector<long double> n_i(domains[col_i], 0);
                std::vector<long double> n_j(domains[col_k], 0);

                ComputeContingencyTable(smp, data, col_i, col_k, n_i_j, n_i, n_j);
                //PrintContingencyTable(n_i_j);
                //FIXME: есть ощущение что можно выставить параметры таким образом, что данная проверка пропустит таблицу
                // для которой вычисление хи - квадрат упадёт (будет деление на 0), возможно надо добавить доп. проверку
                // и в случае деления на 0 считать, что TooMuchStructuralZeroes

                if (TooMuchStructuralZeros(n_i_j, col_i, col_k)) {
                    correlations.emplace_back(col_i, col_k);
                    std::cout << "{" << col_i << "," << col_k << "} Too much structural zeroes" << std::endl
                              << std::endl;
                    continue;
                }

                //Computes chi - squared
                long double chi_squared = ComputeChiSquared(n_i_j, n_i, n_j, col_i, col_k, smp.row_indices.size());

                size_t v = (domains[col_i] - 1) * (domains[col_k] - 1);
                std::cout << "degrees of freedom:" << v << std::endl;
                boost::math::chi_squared dist(v);
                long double t = quantile(dist, 1 - p);
                long double cd = cdf(dist, chi_squared);
                std::cout << "chi squared: " << chi_squared << " quantile: " << t << " cdf: " << cd << std::endl;
                std::cout << std::endl;
                if (chi_squared > t) {
                    correlations.emplace_back(col_i, col_k);
                    continue;
                }
                independent.emplace_back(col_i, col_k);
            }
            //std::cout << "#" << std::endl;
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
        std::cout << "independent: " << std::endl;
        for (auto it: independent) {
            std::cout << "{" << it.first << "," << it.second << "}" << " ";
        }
        std::cout << std::endl;
//        for (auto it: domains) {
//            std::cout << it << " ";
//        }
    }
}