#include <limits>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cassert>
#include "Utils.hpp"
#include "Matrix.hpp"

Mat create_matrix(size_t rows, size_t cols, Allocator *allocator)
{
  Mat res = { (double *)allocate(*allocator,
                                 rows * cols * sizeof (double)),
              rows,
              cols,
              allocator };

  if (res.data == NULL)
  {
    std::fprintf(stderr,
                 "ERROR: failed to allocate matrix with dimensions "
                   "rows x columns = %zu x %zu.\n",
                 rows,
                 cols);
    std::exit(EXIT_FAILURE);
  }

  return res;
}

void destroy(const Mat &matrix)
{
  deallocate(*matrix.allocator, matrix.data);
}

double &at(const Mat &matrix, size_t row, size_t column)
{
  assert(row < matrix.rows && column < matrix.cols);

  return matrix.data[row * matrix.cols + column];
}

Mat copy(const Mat &matrix)
{
  Mat res = create_matrix(matrix.rows,
                          matrix.cols,
                          matrix.allocator);

  std::memcpy(res.data,
              matrix.data,
              matrix.rows * matrix.cols * sizeof (double));

  return res;
}

void fill_randomly(const Mat &matrix)
{
  for (size_t i = matrix.rows * matrix.cols; i-- > 0; )
    matrix.data[i] = rand(-1, 1);
}

Mat multiply(const Mat &left, const Mat &right)
{
  assert(left.cols == right.rows);

  Mat res = create_matrix(left.rows, right.cols, left.allocator);

  for (size_t i = res.rows * res.cols; i-- > 0; )
    res.data[i] = 0;

  for (size_t i = 0; i < left.rows; i++)
  {
    for (size_t j = 0; j < left.cols; j++)
    {
      double const factor = at(left, i, j);
      for (size_t k = 0; k < right.cols; k++)
        at(res, i, k) += factor * at(right, j, k);
    }
  }

  return res;
}

#define PIVOT_CLOSE_TO_ZERO 0x1
#define PIVOT_SWAPPED 0x2

uint8_t pivot(Mat &matrix, size_t column, size_t *perms)
{
  assert(matrix.rows == matrix.cols);

  double best_value = 0;
  size_t best_column = column;

  for (size_t i = column; i < matrix.cols; i++)
  {
    double const val = std::abs(at(matrix, column, i));

    if (val > best_value)
    {
      best_value = val;
      best_column = i;
    }
  }

  uint8_t const is_close_to_zero = PIVOT_CLOSE_TO_ZERO *
    (best_value < std::numeric_limits<double>::epsilon());

  if (column != best_column)
  {
    for (size_t i = 0; i < matrix.rows; i++)
      std::swap(at(matrix, i, column), at(matrix, i, best_column));

    if (perms != NULL)
      std::swap(perms[column], perms[best_column]);

    return is_close_to_zero | PIVOT_SWAPPED;
  }

  return is_close_to_zero;
}

double abs(const Mat &matrix)
{
  assert(matrix.rows == matrix.cols);

  Mat upper = copy(matrix);
  double det = 1;

  for (size_t i = 0; i < upper.cols; i++)
  {
    uint8_t const flags = pivot(upper, i, NULL);

    if (flags & PIVOT_CLOSE_TO_ZERO)
    {
      det = 0;

      break;
    }

    if (flags & PIVOT_SWAPPED)
      det = -det;

    det *= at(upper, i, i);

    for (size_t j = i + 1; j < upper.rows; j++)
    {
      double const factor = -at(upper, j, i) / at(upper, i, i);
      for (size_t k = i + 1; k < upper.cols; k++)
        at(upper, j, k) += factor * at(upper, i, k);
    }
  }

  delete[] upper.data;

  return det;
}

void *malloc_or_panic(size_t size)
{
  void *const data = std::malloc(size);

  if (data == NULL)
    std::exit(EXIT_FAILURE);

  return data;
}

Mat inverse(const Mat &matrix)
{
  assert(matrix.rows == matrix.cols);

  Mat lu = copy(matrix);
  size_t *const perms =
    (size_t *)malloc_or_panic(matrix.cols * sizeof (size_t));

  for (size_t i = 0; i < lu.cols; i++)
    perms[i] = i;

  for (size_t i = 0; i < lu.cols; i++)
  {
    if (pivot(lu, i, perms) & PIVOT_CLOSE_TO_ZERO)
    {
      std::fputs("ERROR: matrix is not invertible.\n", stderr);
      std::exit(EXIT_FAILURE);
    }

    for (size_t j = i + 1; j < lu.rows; j++)
    {
      double const factor = -at(lu, j, i) / at(lu, i, i);

      at(lu, j, i) = factor;
      for (size_t k = i + 1; k < lu.cols; k++)
        at(lu, j, k) += factor * at(lu, i, k);

      // Invert lower matrix.
      for (size_t k = 0; k < i; k++)
        at(lu, j, k) += factor * at(lu, i, k);
    }
  }

  // Invert upper matrix.
  for (size_t i = lu.cols; i-- > 0; )
  {
    for (size_t j = i; j-- > 0; )
    {
      double const factor = -at(lu, j, i) / at(lu, i, i);

      at(lu, j, i) = factor;
      for (size_t k = i + 1; k < lu.cols; k++)
        at(lu, j, k) += factor * at(lu, i, k);
    }

    // Set main diagonal to ones.
    double const factor = 1 / at(lu, i, i);

    at(lu, i, i) = factor;
    for (size_t j = i + 1; j < lu.cols; j++)
      at(lu, i, j) *= factor;
  }

  Mat tmp = { (double *)malloc_or_panic(matrix.rows * matrix.cols * sizeof (double)),
              matrix.rows,
              matrix.cols,
              NULL };

  // Multiply lower and upper matrices.
  for (size_t i = matrix.rows * matrix.cols; i-- > 0; )
    tmp.data[i] = 0;

  for (size_t i = 0; i < lu.rows; i++)
  {
    for (size_t j = i; j < lu.cols; j++)
    {
      double const factor = at(lu, i, j);

      at(tmp, i, j) += factor;
      for (size_t k = 0; k < j; k++)
        at(tmp, i, k) += factor * at(lu, j, k);
    }
  }

  // Permutate product of lower and upper matrices back to its
  // original position.
  for (size_t i = 0; i < matrix.cols; i++)
  {
    std::memcpy(&at(lu, perms[i], 0),
                &at(tmp, i, 0),
                matrix.cols * sizeof (double));
  }

  std::free(tmp.data);
  std::free(perms);

  return lu;
}

void print(const Mat &matrix)
{
  std::putchar('{');

  for (size_t i = 0; i < matrix.rows; i++)
  {
    std::putchar('{');

    for (size_t j = 0; j < matrix.cols; j++)
    {
      std::printf("%f%c",
                  at(matrix, i, j),
                  j + 1 < matrix.cols ? ',' : '}');
    }

    if (i + 1 < matrix.rows)
      std::putchar(',');
  }

  std::putchar('}');
}
