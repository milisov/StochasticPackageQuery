#!/bin/bash

cfg_file="config.cfg"
build_type=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
hostname=$(grep 'hostname' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
port=$(grep 'port' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
username=$(grep 'username' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
database=$(grep 'database' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
password=$(grep 'password' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

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

pip3 install -r requirements.txt
conan profile detect --force
if [ ! -d build ]; then
    mkdir build
fi
conan install . --output-folder=build --build=missing --settings=build_type=$build_type

cd test
echo "Creating table spaqls..."
python3 sqls.py
cd ..

bash test/create_tables.sh