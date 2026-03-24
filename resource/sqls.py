#!/usr/bin/env python
# coding: utf-8

# In[1]:


import configparser, os, glob
spaql_file = '../resource/sqls/{}.spaql'
config_file = '../config.cfg'
config = configparser.ConfigParser()
config.read(config_file)
config['partition'] = {}


# In[2]:


sqls_folder = '../resource/sqls'
pattern = os.path.join(sqls_folder, '_*.spaql')
files_to_remove = glob.glob(pattern)
for file_path in files_to_remove:
    try:
        os.remove(file_path)
    except OSError as e:
        print(f"Error: {file_path} : {e.strerror}")


# In[3]:
stock_tables = []

stock_table = 'stocks'
nstocks = [3, 4, 5, 6]
npaths = [2, 4]
current_stock_tables = [f'{stock_table}_{nstock}_{npath}' for nstock in nstocks for npath in npaths]
with open(spaql_file.format(stock_table), 'r') as file:
    content = file.read()
for table in current_stock_tables:
    modified_content = content.format(table)
    with open(spaql_file.format("_"+table), 'w') as file:
        file.write(modified_content)
    config['partition'][table] = 'price,profit'

stock_tables.extend(current_stock_tables)

stock_table = 'stocks'
nstocks = [3, 4, 5, 6]
npaths = [2, 4]
current_stock_tables = [f'{stock_table}_{nstock}_{npath}_validate' for nstock in nstocks for npath in npaths]
for table in current_stock_tables:
    config['partition'][table] = 'price,profit'

# stock_tables.extend(current_stock_tables)

stock_table = 'stocks'
nstocks = [3,4,5]
npaths = [2, 4]
seeds = list(range(1,11))
current_stock_tables = [f'{stock_table}_{nstock}_{npath}_seeded_{seed}' for nstock in nstocks for npath in npaths for seed in seeds]
with open(spaql_file.format(stock_table), 'r') as file:
    content = file.read()
for table in current_stock_tables:
    modified_content = content.format(table)
    with open(spaql_file.format("_"+table), 'w') as file:
        file.write(modified_content)
    config['partition'][table] = 'price,profit'

stock_tables.extend(current_stock_tables)


stock_table = 'stocks'
nstocks = [3,4,5]
npaths = [2, 4]
seeds = list(range(1,11))
current_stock_tables = [f'{stock_table}_{nstock}_{npath}_seeded_{seed}_validate' for nstock in nstocks for npath in npaths for seed in seeds]
for table in current_stock_tables:
    config['partition'][table] = 'price,profit'


# In[4]:


with open('create_tables.sh', 'w') as bashfile:
    lines = []
    for table in stock_tables:
        lines.append(f'echo "Creating table {table}..."')
        command = ' '.join(table.split('_')[1:])
        lines.append(f"python3 stocks.py {command}")
    bashfile.write('\n'.join(lines))


# In[5]:


with open('../_config.cfg', 'w') as configfile:
    config.write(configfile)

