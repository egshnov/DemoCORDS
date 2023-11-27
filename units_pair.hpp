#pragma once

#include "candidate.h"

namespace algos {
    struct UnitsPair { //TODO:rename
        size_t col_l;
        size_t col_r;
        size_t row_ind;

        UnitsPair(size_t l, size_t r, size_t ind) {
            col_l = l;
            col_r = r;
            row_ind = ind;
        }
    };
}