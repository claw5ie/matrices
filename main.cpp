#include <limits>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cassert>

template<class Type>
Type *alloc(size_t count)
{
  Type *const data = (Type *)std::malloc(count * sizeof (Type));

  if (data == NULL)
  {
    std::fprintf(stderr,
                 "ERROR: failed to allocate %zu elements of "
                   "size %zu.\n",
                 count,
                 sizeof (Type));
    std::exit(EXIT_FAILURE);
  }

  return data;
}

double rand(double min, double max)
{
  return (double)std::rand() / RAND_MAX * (max - min) + min;
}

struct Mat
{
  double *data;
  size_t rows;
  size_t cols;
};

void fill_randomly(Mat &matrix)
{
  for (size_t i = matrix.rows * matrix.cols; i-- > 0; )
    matrix.data[i] = rand(-1, 1);
}

double &at(Mat &matrix, size_t row, size_t column)
{
  assert(row < matrix.rows && column < matrix.cols);

  return matrix.data[row * matrix.cols + column];
}

const double &at(const Mat &matrix, size_t row, size_t column)
{
  assert(row < matrix.rows && column < matrix.cols);

  return matrix.data[row * matrix.cols + column];
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

Mat copy(const Mat &matrix)
{
  Mat res = matrix;

  size_t const count = matrix.rows * matrix.cols;
  res.data = alloc<double>(count);
  std::memcpy(res.data, matrix.data, count * sizeof (double));

  return res;
}

Mat &multiply(Mat &dest, const Mat &left, const Mat &right)
{
  assert(left.cols == right.rows &&
         dest.rows == left.rows &&
         dest.cols == right.cols);

  for (size_t i = dest.rows * dest.cols; i-- > 0; )
    dest.data[i] = 0;

  for (size_t i = 0; i < left.rows; i++)
  {
    for (size_t j = 0; j < left.cols; j++)
    {
      double const factor = at(left, i, j);
      for (size_t k = 0; k < right.cols; k++)
        at(dest, i, k) += factor * at(right, j, k);
    }
  }

  return dest;
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

Mat inverse(const Mat &matrix)
{
  assert(matrix.rows == matrix.cols);

  Mat lu = copy(matrix);
  size_t *const perms = alloc<size_t>(matrix.cols);

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

  Mat tmp = { alloc<double>(matrix.rows * matrix.cols),
              matrix.rows,
              matrix.cols };

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

  delete[] tmp.data;
  delete[] perms;

  return lu;
}

int main()
{
  Mat x = { NULL, 8, 8 };
  x.data = alloc<double>(x.rows * x.cols);
  Mat res = { NULL, x.rows, x.cols };
  res.data = alloc<double>(res.rows * res.cols);

  fill_randomly(x);
  Mat inv = inverse(x);
  print(multiply(res, inv, x));
  std::putchar('\n');

  delete[] x.data;
  delete[] res.data;
  delete[] inv.data;
}
