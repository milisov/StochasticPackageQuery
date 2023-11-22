#!/bin/bash

cfg_file="config.cfg"
BUILD_TYPE=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

# Check if PostgreSQL is installed
if ! command -v pg_config &> /dev/null; then
  echo "PostgreSQL is not installed. Exiting the script."
  exit 1
fi

pip3 install -r requirements.txt
conan profile detect --force
if [ ! -d build ]; then
    mkdir build
fi
conan install . --output-folder=build --build=missing --settings=build_type=$BUILD_TYPE

cd test
echo "Creating table stocks3..."
python3 stocks.py 3
echo "Creating table stocks30..."
python3 stocks.py 30
echo "Creating table stocks300..."
python3 stocks.py 300
echo "Creating table stocks3000..."
python3 stocks.py 3000
cd ..