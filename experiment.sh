#!/bin/bash

. .bashrc

# Default values
N=0
M=0
ALGORITHM=""
H="" # <-- 1. Initialize H with a default (empty string)

# Parse arguments
# 2. Add 'H:' to the getopts string to recognize the -H flag with a value.
while getopts "N:M:a:H:" opt; do
  case $opt in
    N) N=$OPTARG ;;
    M) M=$OPTARG ;;
    a) ALGORITHM=$OPTARG ;;
    H) H=$OPTARG ;; # <-- 3. Add a case to handle the -H flag
    # 4. Update the usage message to include the optional [-H] flag.
    *) echo "Usage: $0 -N <number> -M <number> -a <algorithm> [-H <value>]"; exit 1 ;;
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

# Run the main executable with the provided arguments
# 5. Pass the new "$H" variable as an argument to your main executable.
# If -H was not provided, this will pass an empty string "".
./main "$N" "$M" "$ALGORITHM" "$H"