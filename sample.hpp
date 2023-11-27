#pragma once

#include <vector>
#include "units_pair.hpp"

namespace algos {
    struct Sample {
        std::vector<UnitsPair> self;
        size_t lhs_cardinality;
        size_t rhs_cardinality;
    };
}