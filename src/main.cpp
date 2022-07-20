#include <cstdio>
#include <cstdlib>
#include "Matrix.hpp"

int main()
{
  LinearAllocator linear_allocator = create_linear_allocator(4096);
  Allocator allocator;

  allocator.type = Allocator::Linear;
  allocator.as.linear_allocator = &linear_allocator;

  Mat x = create_matrix(8, 8, &allocator);

  fill_randomly(x);

  print(multiply(inverse(x), x));
  std::putchar('\n');

  std::free(linear_allocator.data);
}
