#pragma once

#include <vector>
#include <string>
#include <fstream>
#include "sample.hpp"

namespace algos {

    class Cords {
    public:
    private:
        bool has_header = false;

        const double pi = 3.14159265358979323846; //TODO: CHANGE to math.pi
        //algo params:
        long double eps1 = 0.04; // (1-eps1) * |R| (where |R| is amount of rows in table R) denotes minimum cardinality of a column to be considered a soft key
        long double eps2 = 0.2; // pre-check for detection of SFD. For deep explanation look at "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies." step in Part "Detection algorithm"
        long double eps3 = 0.3; // (1-eps3) represents minimal strength of a soft FD in a sample
        long double eps4 = 0.1; // (1-eps4) * |R| represents minimal frequency of a value that states that value distribution in a column is skewed.
        long double eps5 = 1e-01; // eps5 * domain1 * domain2 represents minimum amount of structural zeroes in a contingency table that states that columns are correlated
        long double p = 1e-06; // max allowed probability of a false positive result for the chi - squared test
        long double delta = 0.05; // fixed constant. Note that delta must be more than eps1. For deep explanation look at "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies." page 6
        size_t num_freqvals = 70; //represents amount of values which frequencies are tracked to check whether data distribution is skewed and are used in chi-squared test
         long long d_max = 70; //max amount of categories for the chi-squared test in case the data is not skewed

        std::vector<long long> cardinality;
        std::vector<std::unordered_map<std::string, std::pair<size_t, size_t>>> frequent_values_statistics; //pair.first = frequency pair.second = frequency sequence number
        std::vector<size_t> freq_sums;
        std::vector<bool> is_skewed;
        std::vector<long long> domains;

        //Results:
        std::vector<size_t> soft_keys{};
        std::vector<size_t> trivial_columns{};
        std::vector<std::pair<size_t, size_t>> correlations;
        std::vector<std::pair<size_t, size_t>> independent;
        std::vector<std::pair<size_t, size_t>> sfd; // first represents C1, second represents C2 where C1->C2

        void GetFrequentValuesStatistics(const Table &data);

        void
        Filter(const std::unordered_map<std::string, std::pair<size_t, size_t>> &freq_val_stats, const Table &data,
               size_t col_ind,
               Sample &smp);

        //TODO: Determine val type before adding to desbordante (most likely Type*)
        [[nodiscard]] size_t
        Category(size_t col_ind, const std::string &val, size_t domain, bool skew, const Table &data) const;

        //TODO: separate struct for contingency table
        void ComputeContingencyTable(Sample &smp, const Table &data, size_t col_i, size_t col_k,
                                     std::vector<std::vector<long double>> &n_i_j, std::vector<long double> &n_i,
                                     std::vector<long double> &n_j);

        bool TooMuchStructuralZeros(const std::vector<std::vector<long double>> &n_i_j, size_t col_i, size_t col_k);

        long double
        ComputeChiSquared(const std::vector<std::vector<long double>> &n_i_j, const std::vector<long double> &n_i,
                          const std::vector<long double> &n_j, size_t col_i, size_t col_k, long double sample_size);

        void Init(size_t columns);

        std::pair<size_t, size_t> SwitchIndicesIfNecessary(size_t ind1, size_t ind2);


        std::vector<size_t> GetSoftKeys() {
            return soft_keys;
        }

        [[nodiscard]] long long ComputeSampleSize(long long d1, long long d2) const;

        //Sample Sampling(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table &data);
        bool DetectSfd(const Sample &smp);

        void SkewHandling(size_t col_i, size_t col_k, const Table &data, Sample &smp);

        static Table GetTable(const std::string &filename);

    public:
        Cords();

        void ExecuteInternal(const std::string &filename, bool has_header);
    };
}