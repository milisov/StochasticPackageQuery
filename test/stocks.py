#!/usr/bin/env python
# coding: utf-8

# In[16]:


import QuantLib as ql
import numpy as np
import matplotlib.pyplot as plt
import configparser
import multiprocessing
import psycopg2
import os, random, sys
from psycopg2 import Error


# In[17]:


nStocks = 3000
nPaths = 10000


# In[ ]:


if __name__ == "__main__":
    if len(sys.argv) > 1:
        nStocks = int(sys.argv[1])
    if len(sys.argv) > 2:
        nPaths = int(sys.argv[2])


# In[15]:


config = configparser.ConfigParser()
config.read('../config.cfg')
is_rebuild = config['build']['rebuild_stocks'] == 'true'
maturity = 1.0
nSteps = int(maturity * 365)
timeGrid = np.linspace(0.0, maturity, nSteps + 1)

pairs = []
num_cores = multiprocessing.cpu_count() // 2
table_name = "stocks"+str(nStocks)

stocks = configparser.ConfigParser()
stocks.read('../resource/stocks/tickers.ini')
stats = []
for ticker in stocks.sections():
    stats.append([ticker, float(stocks[ticker]['price']), float(stocks[ticker]['volatility']), float(stocks[ticker]['drift'])])


# In[3]:


def GeneratePaths(process, maturity, nPaths, nSteps):
    generator = ql.UniformRandomGenerator(42)
    sequenceGenerator = ql.UniformRandomSequenceGenerator(nSteps, generator)
    gaussianSequenceGenerator = ql.GaussianRandomSequenceGenerator(sequenceGenerator)
    paths = np.zeros(shape = (nPaths, nSteps+1))
    pathGenerator = ql.GaussianPathGenerator(process, maturity, nSteps, gaussianSequenceGenerator, False)
    for i in range(nPaths):
        path = pathGenerator.next().value()
        paths[i, :] = np.array(path)
    return paths

# Define your procedure here
def simulate(stat):
    ticker, price, vol, drift = stat
    GBM = ql.GeometricBrownianMotionProcess(price, vol, drift)
    gbm_paths = GeneratePaths(GBM, maturity, nPaths, nSteps)[:, 1:] - price
    process = multiprocessing.current_process()
    process_id = (process._identity[0]-1)%num_cores
    conn, cur = pairs[process_id]
    for day_index in range(nSteps):
        try:
            cur.execute(f"""
                INSERT INTO {table_name} (stock, price, day, profit)
                VALUES (%s, %s, %s, %s)
            """, (ticker, price, day_index+1, list(gbm_paths[:, day_index])))
        except psycopg2.Error as e:
            continue
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


# In[4]:


conn, cur = get_conn_cur(config)
if table_exists(cur, table_name):
    if is_rebuild:
        delete_table_sql = f"DROP TABLE IF EXISTS {table_name}"
        cur.execute(delete_table_sql)
        conn.commit()
is_populating = False
if not table_exists(cur, table_name):
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id SERIAL PRIMARY KEY,
            stock TEXT,
            price REAL,
            day INT,
            profit REAL[]
        )
    """)
    conn.commit()
    is_populating = True
cur.close()
conn.close()


# In[5]:


if is_populating:
    for i in range(num_cores):
        pairs.append(get_conn_cur(config))
    pool = multiprocessing.Pool(processes=num_cores)
    random.seed(42)
    random.shuffle(stats)
    pool.map(simulate, stats[:nStocks])
    pool.close()
    pool.join()
    for conn, cur in pairs:    
        cur.close()
        conn.close()


# In[ ]:




