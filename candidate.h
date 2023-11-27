#pragma once

#include <memory>

namespace algos {

    struct Candidate {
        /* Column indexes, the first for the column whose values were the left operand
         * for binop_, the second for the right */
        size_t l;
        size_t r;
    };

}  // namespace algos
