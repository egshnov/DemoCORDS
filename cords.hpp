#pragma once

#include <vector>
#include <string>
#include <fstream>
#include "sample.hpp"

namespace algos {
    class Cords {
    public:
        using Table = std::vector<std::vector<std::string>>;
    private:
        const double pi = 3.14159265358979323846;
        long double eps1 = 0.04; // (1-eps1) * |R| where |R| is amount of rows in table R denotes minimum cardinality of a column to be considered a soft key
        long double eps2 = 0.3;
        long double eps3 = 0.2;
        long double p = 1e-06;// max allowed probability of a false positive result for the chi - squared test
        long double delta = 0.05; // fixed constant. Note that delta must be more than eps1. For deep explanation look
        // at "CORDS: Automatic Discovery of Correlations and Soft Functional Dependencies." page 6
        std::vector<size_t> soft_keys{};

        unsigned long long ComputeCardinality(size_t ind, const Table &data);

        //unsigned long long ComputeSampleColumnCardinality(size_t ind, const Table &data, const std::vector<UnitsPair> &sample);

    public:
        Cords();

        std::vector<size_t> GetSoftKeys() {
            return soft_keys;
        }

        unsigned long long ComputeSampleSize(unsigned long long d1, unsigned long long d2) const;
        //TODO: сделать конструктором для Sample ???
        Sample Sampling(unsigned long long sample_size, size_t rows, size_t l, size_t r, const Table& data);

        static Table GetTable(const std::string &filename);

        void ExecuteInternal(const std::string &filename);
    };
}