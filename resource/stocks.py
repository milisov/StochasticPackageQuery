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
        except ValueError:
            pass


# In[4]:


nStocks = 10**nStocksLog
nPaths = 10**nPathsLog


# In[5]:


config = configparser.ConfigParser()
config.read('../config.cfg')
is_rebuild = config['build']['rebuild_stocks'] == 'true'
validate_scenarios = int(config['build']['validate_scenarios'])
validate_seed = int(config['build']['validate_seed'])
generate_seed = int(config['build']['generate_seed'])
maturity = 1.0
stepPerYear = 52
nSteps = int(maturity * stepPerYear)

pairs = []
num_cores = multiprocessing.cpu_count() // 2
table_name = f"stocks_{nStocksLog}_{nPathsLog}"
validate_table_name = f"stocks_{nStocksLog}_{nPathsLog}_validate"

stocks = configparser.ConfigParser()
stocks.read('../resource/stocks/tickers.ini')
stats = []
for ticker in stocks.sections():
    stats.append([ticker, float(stocks[ticker]['price']), float(stocks[ticker]['volatility']), float(stocks[ticker]['drift'])])
POS_INF = 1e30
NEG_INF = -1e30


# In[6]:


def GeneratePaths(process, maturity, nPaths, nSteps, isValidate):
    if isValidate:
        generator = ql.UniformRandomGenerator(validate_seed)
    else:
        generator = ql.UniformRandomGenerator(generate_seed)
    sequenceGenerator = ql.UniformRandomSequenceGenerator(nSteps, generator)
    gaussianSequenceGenerator = ql.GaussianRandomSequenceGenerator(sequenceGenerator)
    paths = np.zeros(shape = (nPaths, nSteps+1), dtype=np.float32)
    pathGenerator = ql.GaussianPathGenerator(process, maturity, nSteps, gaussianSequenceGenerator, False)
    for i in range(nPaths):
        path = pathGenerator.next().value()
        paths[i, :] = np.clip(np.array(path), NEG_INF, POS_INF)
    return paths

# Define your procedure here
def simulate(stat):
    start_id, (ticker, price, vol, drift) = stat
    process = multiprocessing.current_process()
    process_id = (process._identity[0]-1)%num_cores
    conn, cur = pairs[process_id]
    GBM = ql.GeometricBrownianMotionProcess(price, vol, drift)
    gbm_paths = GeneratePaths(GBM, maturity, nPaths, nSteps, False)[:, 1:] - price
    data = io.StringIO()
    for week_index in range(nSteps):
        profit = "{" + ",".join(map(str, gbm_paths[:, week_index])) + "}"
        day_index = (week_index+1)*(365//stepPerYear)
        data.write(f"{start_id+week_index}|'{ticker}'|{day_index}|{price}|{profit}\n")
    data.seek(0)
    cur.copy_from(data, table_name, sep='|')
    validate_gbm_paths = GeneratePaths(GBM, maturity, validate_scenarios, nSteps, True)[:, 1:] - price
    validate_data = io.StringIO()
    for week_index in range(nSteps):
        profit = "{" + ",".join(map(str, validate_gbm_paths[:, week_index])) + "}"
        day_index = (week_index+1)*(365//stepPerYear)
        validate_data.write(f"{start_id+week_index}|'{ticker}'|{day_index}|{price}|{profit}\n")
    validate_data.seek(0)
    cur.copy_from(validate_data, validate_table_name, sep='|')
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
is_populating = False
def check_table(cur, table):
    if table_exists(cur, table):
        if is_rebuild:
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
is_populating = check_table(cur, table_name)
is_populating = check_table(cur, validate_table_name)
cur.close()
conn.close()


# In[8]:


if is_populating:
    for i in range(num_cores):
        pairs.append(get_conn_cur(config))
    pool = multiprocessing.Pool(processes=num_cores)
    random.seed(generate_seed)
    random.shuffle(stats)
    nStocks = min(nStocks//stepPerYear, len(stats))
    subStats = sorted(stats[:nStocks])
    start_ids = [i*stepPerYear+1 for i in range(nStocks)]
    pool.map(simulate, list(zip(start_ids, subStats)))
    pool.close()
    pool.join()
    for conn, cur in pairs:    
        cur.close()
        conn.close()
    conn, cur = get_conn_cur(config)
    cur.execute(f"CREATE INDEX id_{table_name} ON {table_name} (id);")
    cur.execute(f"CREATE INDEX id_{validate_table_name} ON {validate_table_name} (id);")
    conn.commit()
    cur.close()
    conn.close()


# In[ ]:




