#!/bin/bash

cfg_file="config.cfg"
gb_key=$(grep 'gurobi_key' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
build_type=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
hostname=$(grep 'hostname' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
port=$(grep 'port' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
username=$(grep 'username' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
database=$(grep 'database' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
password=$(grep 'password' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

project_folder=$(pwd)

if [ ! -f "$project_folder/.bashrc" ]; then
    touch "$project_folder/.bashrc"
fi

if ! command -v cmake >/dev/null 2>&1; then
    if [ ! -d resource/cmake ]; then
        mkdir resource/cmake
    fi
    cd resource/cmake
    cmake_name="cmake-3.28.0-rc5-linux-x86_64"
    if [ ! -f "$cmake_name.sh" ]; then
        wget https://github.com/Kitware/CMake/releases/download/v3.28.0-rc5/$cmake_name.sh
    fi
    if [ ! -d "$cmake_name" ]; then
        chmod +x $cmake_name.sh
        echo -e "y\ny" | ./$cmake_name.sh --prefix=.
    fi
    cd $cmake_name
    if ! grep -q "$cmake_name" $project_folder/.bashrc; then
        echo "export PATH=$(pwd)/bin:\$PATH" >> $project_folder/.bashrc
    fi
    cd ../../..
else
    cmp=3.21
    ver=$(cmake --version | head -1 | cut -f3 -d" ")
    mapfile -t sorted < <(printf "%s\n" "$ver" "$cmp" | sort -V)
    if ! [[ ${sorted[0]} == "$cmp" ]]; then
        echo "cmake version needs to be at least $cmp"
        exit 1
    fi
fi

if [ ! -d resource/gurobi ]; then
    mkdir resource/gurobi
fi
gb_install_file="gurobi11.0.0_linux64.tar.gz"
gb_version="1100"
cd resource/gurobi
if [ ! -f "$gb_install_file" ]; then
    curl -O https://packages.gurobi.com/11.0/$gb_install_file
fi
if [ ! -d "gurobi$gb_version" ]; then
    tar xvfz $gb_install_file
fi
cd gurobi$gb_version/linux64
gb_install_dir=$(pwd)
if ! grep -q "GUROBI_HOME" $project_folder/.bashrc; then
    echo "export GUROBI_HOME=$gb_install_dir" >> $project_folder/.bashrc
    echo "export PATH=$gb_install_dir/bin:\$PATH" >> $project_folder/.bashrc
    echo "export LD_LIBRARY_PATH=$gb_install_dir/lib:\$LD_LIBRARY_PATH" >> $project_folder/.bashrc
fi
cd ../..
if [ ! -f "gurobi.lic" ]; then
    echo . | $gb_install_dir/bin/grbgetkey $gb_key
fi
if ! grep -q "GRB_LICENSE_FILE" $project_folder/.bashrc; then
    echo "export GRB_LICENSE_FILE=$(pwd)/gurobi.lic" >> $project_folder/.bashrc
fi
cd ../..

# Check if psql command is available
if ! command -v psql &> /dev/null
then
    echo "PostgreSQL(psql) is not installed."
    exit 1
fi

# Get the version of PostgreSQL
version=$(psql --version | awk '{print $3}')
major_version=$(echo $version | cut -d'.' -f1)

# Check if the version is greater than or equal to 14
if [ "$major_version" -ge 14 ]; then
    echo "PostgreSQL version $version is installed."
else
    echo "PostgreSQL version $version is installed, but it is less than 14."
    exit 1
fi

total_cores=$(nproc)
echo "Number of logical cores:" $total_cores

export PGPASSWORD="$password"

if [ -z "$username" ]; then
    max_connections=$(psql -h $hostname -p $port -d $database -t -c "SHOW max_connections;")
else
    max_connections=$(psql -h $hostname -p $port -U $username -d $database -t -c "SHOW max_connections;")
fi

if [ $((total_cores*2)) -ge $max_connections ]; then
    echo "Requires PostgreSQL's max_connections to be least twice the number of logical cores."
    echo "You can change PostgreSQL's max_connections via '/etc/postgresql/<version>/main/postgresql.conf'."
    echo "Remember to change PostgreSQL's shared_buffers accordingly to handle the increase in max_connections."
    echo "And then restart the server using 'sudo systemctl restart postgresql'."
    exit 1
fi

#pip3 install -r requirements.txt
# conan profile detect --force
if [ -d build ]; thens
    rm -r build
fi
conan install . --output-folder=. --build=missing --settings=build_type=$build_type

# echo "Creating table spaqls..."
# cd resource
# python3 sqls.py
# bash create_tables.sh
# cd ..