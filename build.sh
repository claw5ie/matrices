#!/bin/bash

set -xeu

warning_flags="-Wall -Wextra -pedantic"
other_flags="-g -std=c++11"
files="src/main.cpp src/Allocators.cpp src/Matrix.cpp src/Utils.cpp"

g++ ${warning_flags} ${other_flags} ${files}
