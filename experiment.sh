#!/bin/bash

. .bashrc

# Default values
N=0
M=0
ALGORITHM=""

# Parse arguments
while getopts "N:M:a:" opt; do
  case $opt in
    N) N=$OPTARG ;;
    M) M=$OPTARG ;;
    a) ALGORITHM=$OPTARG ;;
    *) echo "Usage: $0 -N <number> -M <number> -a <algorithm>"; exit 1 ;;
  esac
done

cfg_file="config.cfg"
BUILD_TYPE=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
REBUILD=$(grep 'rebuild_cmake' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

cd resource
python3 sqls.py
cd ..

cd build/$BUILD_TYPE/generators
if [ -f "CMakeCache.txt" ]; then
  [ "$REBUILD" = true ] && rm CMakeCache.txt && cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ../../..
else
  cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ../../..
fi
cmake --build . --config $BUILD_TYPE

# Run the main executable with the provided arguments  //make a statement if !RCL --> c++, else python script with respective flag
./main "$N" "$M" "$ALGORITHM"
