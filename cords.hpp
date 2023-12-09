#pragma once

#include <vector>
#include <string>
#include <fstream>
#include "sample.hpp"
#include "alliacies.hpp"

namespace algos {
    class Cords {
    public:
    private:
        bool has_header = false;

        const double pi = 3.14159265358979323846;
        long double eps1 = 0.04; // (1-eps1) * |R| where |R| is amount of rows in table R denotes minimum cardinality of a column to be considered a soft key
        long double eps2 = 0.2; // (1-eps2) represents required proximity of a |C1C2|_S to |S|
        long double eps3 = 0.3; // (1-eps3) represents minimal strength of a soft FD in a sample
        long double eps4 = 0.1;// (1-eps4) * |R| represents
        long double eps5 = 1e-06;
        long double p = 1e-06;// max allowed probability of a false positive result for the chi - squared test
        long double delta = 0.05; // fixed constant. Note that delta must be more than eps1. For deep explanation look at "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies." page 6
        size_t num_freqvals = 70; //represents amount of values which frequencies are tracked to check whether data distribution is skewed
        unsigned long long d_max = 50; //max amount of categories for the chi-squared test in case the data is not skewed

        std::vector<unsigned long long> cardinality;
        std::vector<std::unordered_map<std::string, std::pair<size_t, size_t>>> frequent_values_statistics; //pair.first = frequency pair.second = frequency sequence number
        std::vector<size_t> freq_sums;
        std::vector<bool> is_skewed;
        std::vector<size_t> domains;


        std::vector<size_t> soft_keys{};
        std::vector<size_t> trivial_columns{};
        std::vector<std::pair<size_t, size_t>> correlations;
        std::vector<std::pair<size_t, size_t>> sfd; // first represents C1 second represents C2 where C1->C2

        void GetFrequentValuesStatistics(const Table &data);

        void
        Filter(const std::unordered_map<std::string, std::pair<size_t, size_t>> &freq_val_stats, const Table &data,
               size_t col_ind,
               Sample &smp);

        //TODO: Determine val type
        [[nodiscard]] size_t
        Category(size_t col_ind, const std::string &val, size_t domain, bool skew, const Table &data) const;

    public:
        Cords();

        std::vector<size_t> GetSoftKeys() {
            return soft_keys;
        }

        [[nodiscard]] long long ComputeSampleSize(unsigned long long d1, unsigned long long d2) const;

        //TODO: сделать конструктором для Sample ???
        //Sample Sampling(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table &data);
        bool DetectSfd(const Sample &smp);

        static Table GetTable(const std::string &filename);

        void ExecuteInternal(const std::string &filename, bool has_header);
    };
}