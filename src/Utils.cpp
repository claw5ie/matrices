#include <cstdlib>
#include "Utils.hpp"

double rand(double min, double max)
{
  return (double)std::rand() / RAND_MAX * (max - min) + min;
}
