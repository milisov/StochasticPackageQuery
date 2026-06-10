#!/usr/bin/env python
# coding: utf-8

# In[1]:


import QuantLib as ql
import numpy as np
import matplotlib.pyplot as plt
import configparser
import multiprocessing
import psycopg2
import os, random, sys, io
from psycopg2 import Error


# In[2]:


nStocksLog = 3
nPathsLog = 2
seeded = False
seedVal = None

# In[3]:


if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            nStocksLog = int(sys.argv[1])
        except ValueError:
            pass
    if len(sys.argv) > 2:
        try:
            nPathsLog = int(sys.argv[2])
            if len(sys.argv) > 4:
                seeded = (sys.argv[3] == 'seeded')
                seedVal = int(sys.argv[4])
                print(f"I am setting a seed value of {seedVal}")
        except ValueError:
            pass


# In[4]:


nStocks = 10**nStocksLog
nPaths = nPathsLog  # nPathsLog is now the actual number of scenarios, not a power of 10


# In[5]:


config = configparser.ConfigParser()
config.read('../config.cfg')
is_rebuild = config['build']['rebuild_stocks'] == 'true'
validate_scenarios = int(config['build']['validate_scenarios'])
validate_seed = int(config['build']['validate_seed'])
generate_seed = int(config['build']['generate_seed'])
partition_seed = int(config['build']['partition_seed'])
stochastic_seed = int(config['build']['generate_seed'])
if seeded:
    stochastic_seed = seedVal
maturity = 1.0
stepPerYear = 52
nSteps = int(maturity * stepPerYear)

pairs = []
num_cores = multiprocessing.cpu_count() - 10
print("CORES:",num_cores)
table_name = None
validate_table_name = None
partition_table_name = None
if not seeded:
    table_name = f"stocks_{nStocksLog}_{nPathsLog}"
    validate_table_name = f"stocks_{nStocksLog}_validate"
    partition_table_name = f"stocks_{nStocksLog}_partition"
else:
    table_name = f"stocks_{nStocksLog}_{nPathsLog}_seeded_{seedVal}"
    validate_table_name = f"stocks_{nStocksLog}_validate"
    partition_table_name = f"stocks_{nStocksLog}_partition"

stocks = configparser.ConfigParser()
stocks.read('../resource/stocks/tickers.ini')
stats = []
for ticker in stocks.sections():
    stats.append([ticker, float(stocks[ticker]['price']), float(stocks[ticker]['volatility']), float(stocks[ticker]['drift'])])
POS_INF = 1e30
NEG_INF = -1e30


# In[6]:


def GeneratePaths(process, maturity, nPaths, nSteps, seed_type):
    if seed_type == 'validate':
        generator = ql.UniformRandomGenerator(validate_seed)
    elif seed_type == 'partition':
        generator = ql.UniformRandomGenerator(partition_seed)
    else:
        generator = ql.UniformRandomGenerator(stochastic_seed) ## Filip changed this from generate_seed to stochastic_seed --> if we use different seeds we want fixed stocks but different scenarios
    sequenceGenerator = ql.UniformRandomSequenceGenerator(nSteps, generator)
    gaussianSequenceGenerator = ql.GaussianRandomSequenceGenerator(sequenceGenerator)
    paths = np.zeros(shape = (nPaths, nSteps+1), dtype=np.float32)
    pathGenerator = ql.GaussianPathGenerator(process, maturity, nSteps, gaussianSequenceGenerator, False)
    for i in range(nPaths):
        path = pathGenerator.next().value()
        paths[i, :] = np.clip(np.array(path), NEG_INF, POS_INF)
    return paths

# Define your procedure here
def simulate(start_id, stat, is_populating_main, is_populating_validate, is_populating_partition):
    ticker, price, vol, drift = stat
    process = multiprocessing.current_process()
    process_id = (process._identity[0]-1)%num_cores
    conn, cur = pairs[process_id]
    GBM = ql.GeometricBrownianMotionProcess(price, vol, drift)

    if is_populating_main:
        # print("Generating Optimization Table")
        gbm_paths = GeneratePaths(GBM, maturity, nPaths, nSteps, 'main')[:, 1:] - price
        data = io.StringIO()
        for week_index in range(nSteps):
            profit = "{" + ",".join(map(str, gbm_paths[:, week_index])) + "}"
            day_index = (week_index+1)*(365//stepPerYear)
            data.write(f"{start_id+week_index}|'{ticker}'|{day_index}|{price}|{profit}\n")
        data.seek(0)
        cur.copy_from(data, table_name, sep='|')

    if is_populating_validate:
        # print("Generating Validation Table")
        validate_gbm_paths = GeneratePaths(GBM, maturity, validate_scenarios, nSteps, 'validate')[:, 1:] - price
        validate_data = io.StringIO()
        for week_index in range(nSteps):
            profit = "{" + ",".join(map(str, validate_gbm_paths[:, week_index])) + "}"
            day_index = (week_index+1)*(365//stepPerYear)
            validate_data.write(f"{start_id+week_index}|'{ticker}'|{day_index}|{price}|{profit}\n")
        validate_data.seek(0)
        cur.copy_from(validate_data, validate_table_name, sep='|')

    if is_populating_partition:
        # print("Generating Partition Table")
        partition_gbm_paths = GeneratePaths(GBM, maturity, validate_scenarios, nSteps, 'partition')[:, 1:] - price
        partition_data = io.StringIO()
        for week_index in range(nSteps):
            profit = "{" + ",".join(map(str, partition_gbm_paths[:, week_index])) + "}"
            day_index = (week_index+1)*(365//stepPerYear)
            partition_data.write(f"{start_id+week_index}|'{ticker}'|{day_index}|{price}|{profit}\n")
        partition_data.seek(0)
        cur.copy_from(partition_data, partition_table_name, sep='|')

    conn.commit()

def table_exists(cur, table_name):
    cur.execute("SELECT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = %s)", (table_name,))
    return cur.fetchone()[0]

def get_conn_cur(config):
    # Connect to the PostgreSQL database
    conn = psycopg2.connect(
        host=config['postgres']['hostname'],
        port=config['postgres']['port'],
        database=config['postgres']['database'],
        user=config['postgres']['username'],
        password=config['postgres']['password']
    )
    # Create a cursor
    cur = conn.cursor()
    return conn, cur


# In[7]:


conn, cur = get_conn_cur(config)
# is_populating = False
def check_table(cur, table):
    if table_exists(cur, table):
        print("Table", table, "already exists.")
        if is_rebuild:
            print("REBUILDING")
            delete_table_sql = f"DROP TABLE IF EXISTS {table}"
            cur.execute(delete_table_sql)
            conn.commit()
    if not table_exists(cur, table):
        cur.execute(f"""
            CREATE TABLE IF NOT EXISTS \"{table}\" (
                id BIGINT,
                stock TEXT,
                day INT,
                price REAL,
                profit REAL[]
            )
        """)
        conn.commit()
        return True
    return False
is_populating_main = check_table(cur, table_name)
is_populating_validate = check_table(cur, validate_table_name)
is_populating_partition = check_table(cur, partition_table_name)
cur.close()
conn.close()


# In[8]:


if is_populating_main or is_populating_validate or is_populating_partition:
    print("populating")
    for i in range(num_cores):
        pairs.append(get_conn_cur(config))
    pool = multiprocessing.Pool(processes=num_cores)
    random.seed(generate_seed)
    random.shuffle(stats)
    nStocks = min(nStocks//stepPerYear, len(stats))
    subStats = sorted(stats[:nStocks])
    start_ids = [i*stepPerYear+1 for i in range(nStocks)]

    # Bundle the start_id, stat, and our flags together for the worker
    tasks = [(s_id, stat, is_populating_main, is_populating_validate, is_populating_partition) for s_id, stat in zip(start_ids, subStats)]

    pool.starmap(simulate, tasks)
    #pool.map(simulate, list(zip(start_ids, subStats)))
    pool.close()
    pool.join()
    for conn, cur in pairs:
        cur.close()
        conn.close()
    conn, cur = get_conn_cur(config)
    if is_populating_main:
        cur.execute(f"CREATE INDEX id_{table_name} ON {table_name} (id);")
    if is_populating_validate:
        cur.execute(f"CREATE INDEX id_{validate_table_name} ON {validate_table_name} (id);")
    if is_populating_partition:
        cur.execute(f"CREATE INDEX id_{partition_table_name} ON {partition_table_name} (id);")
    conn.commit()
    cur.close()
    conn.close()


# In[ ]:




