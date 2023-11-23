#!/usr/bin/env python
# coding: utf-8

# In[30]:


import configparser
spaql_file = '../resource/sqls/{}.spaql'
config_file = '../config.cfg'
config = configparser.ConfigParser()
config.read(config_file)
config['partition'] = {}


# In[31]:


stock_table = 'stocks'
nstocks = [3, 30, 300, 3000]
npaths = [100, 10000]
stock_tables = [f'{stock_table}_{nstock}_{npath}' for nstock in nstocks for npath in npaths]
with open(spaql_file.format(stock_table), 'r') as file:
    content = file.read()
for table in stock_tables:
    modified_content = content.format(table)
    with open(spaql_file.format("_"+table), 'w') as file:
        file.write(modified_content)
    config['partition'][table] = 'price,profit'


# In[35]:


with open('create_tables.sh', 'w') as bashfile:
    lines = ["cd test"]
    for table in stock_tables:
        lines.append(f'echo "Creating table {table}..."')
        command = ' '.join(table.split('_')[1:])
        lines.append(f"python3 stocks.py {command}")
    lines.append("cd ..")
    bashfile.write('\n'.join(lines))


# In[36]:


with open('../_config.cfg', 'w') as configfile:
    config.write(configfile)


# In[ ]:




