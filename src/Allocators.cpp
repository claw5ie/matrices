#include <cstdio>
#include <cstdlib>
#include "Allocators.hpp"

LinearAllocator create_linear_allocator(size_t capacity)
{
  void *const data = std::malloc(capacity);

  if (data == NULL)
  {
    std::fprintf(stderr,
                 "ERROR: failed to allocate %zu bytes for linear "
                   "allocator.\n",
                 capacity);
    std::exit(EXIT_FAILURE);
  }

  return { (char *)data, 0, capacity };
}

void *allocate(LinearAllocator &allocator, size_t size)
{
  // Align "size" by the size of "void *".
  size += sizeof (void *) - 1;
  size -= size % sizeof (void *);

  if (allocator.size + size > allocator.capacity)
    return NULL;

  void *const data = (void *)(allocator.data + allocator.size);
  allocator.size += size;

  return data;
}

void *allocate(Allocator &allocator, size_t size)
{
  void *data = NULL;

  switch (allocator.type)
  {
  case Allocator::Standard:
    data = std::malloc(size);
    break;
  case Allocator::Linear:
    data = allocate(*allocator.as.linear_allocator, size);
    break;
  }

  return data;
}

void deallocate(Allocator &allocator, void *data)
{
  switch (allocator.type)
  {
  case Allocator::Standard:
    std::free(data);
    break;
  case Allocator::Linear:
    break;
  }
}
