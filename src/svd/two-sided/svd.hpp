#pragma once

#include "../../utils/types.hpp"

size_t svd_blocked(struct matrix_t X, struct matrix_t B, struct matrix_t U, struct matrix_t V, size_t block_size);
