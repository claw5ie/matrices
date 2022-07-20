#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include "Allocators.hpp"

struct Mat
{
  double *data;
  size_t rows;
  size_t cols;
  Allocator *allocator;
};

Mat create_matrix(size_t rows, size_t cols, Allocator *allocator);

void destroy(const Mat &matrix);

double &at(const Mat &matrix, size_t row, size_t column);

Mat copy(const Mat &matrix);

void fill_randomly(const Mat &matrix);

Mat multiply(const Mat &left, const Mat &right);

double abs(const Mat &matrix);

Mat inverse(const Mat &matrix);

void print(const Mat &matrix);

#endif // MATRIX_HPP
