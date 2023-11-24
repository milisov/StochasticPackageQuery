#!/bin/bash

cfg_file="config.cfg"
BUILD_TYPE=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
REBUILD=$(grep 'rebuild_cmake' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

cd test
python3 sqls.py
cd ..

cd build
if [ -f "CMakeCache.txt" ]; then
  [ "$REBUILD" = true ] && rm CMakeCache.txt
fi
[ "$REBUILD" = true ] && cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ..
cmake --build . --config $BUILD_TYPE

./main