#ifndef ALLOCATORS_HPP
#define ALLOCATORS_HPP

#include <cstddef>

struct LinearAllocator
{
  char *data;
  size_t size;
  size_t capacity;
};

LinearAllocator create_linear_allocator(size_t capacity);

void *allocate(LinearAllocator &allocator, size_t size);

struct Allocator
{
  enum Type
  {
    Standard,
    Linear
  };

  Type type;

  union
  {
    LinearAllocator *linear_allocator;
  } as;
};

void *allocate(Allocator &allocator, size_t size);

void deallocate(Allocator &allocator, void *data);

#endif // ALLOCATORS_HPP
