#!/usr/bin/env python
# coding: utf-8

import configparser
import psycopg2
import sys

config = configparser.ConfigParser()
config.read('../config.cfg')

def get_conn_cur(config):
    conn = psycopg2.connect(
        host=config['postgres']['hostname'],
        port=config['postgres']['port'],
        database=config['postgres']['database'],
        user=config['postgres']['username'],
        password=config['postgres']['password']
    )
    cur = conn.cursor()
    return conn, cur

def get_violating_ids(cur, validate_table):
    """Get IDs where std(profit) > 6500 from the validate table"""
    cur.execute(f"""
        SELECT t.id
        FROM {validate_table} t, unnest(t.profit) AS x
        GROUP BY t.id
        HAVING stddev(x) > 6500
    """)
    return [row[0] for row in cur.fetchall()]

def get_related_tables(cur, prefix):
    """Get all tables matching stocks_X_* pattern"""
    cur.execute("""
        SELECT table_name FROM information_schema.tables
        WHERE table_schema = 'public' AND table_name LIKE %s
    """, (f"{prefix}%",))
    return [row[0] for row in cur.fetchall()]

def filter_table(conn, cur, table, violating_ids):
    """Remove violating IDs from a table"""
    if not violating_ids:
        print(f"[FILTER] {table}: No violating IDs to remove")
        return

    cur.execute(f"SELECT COUNT(*) FROM {table}")
    count_before = cur.fetchone()[0]

    # Delete rows with violating IDs
    cur.execute(f"""
        DELETE FROM {table}
        WHERE id = ANY(%s)
    """, (violating_ids,))
    conn.commit()

    cur.execute(f"SELECT COUNT(*) FROM {table}")
    count_after = cur.fetchone()[0]
    removed = count_before - count_after
    print(f"[FILTER] {table}: Removed {removed} rows (before: {count_before}, after: {count_after})")

def reindex_table(conn, cur, table):
    """Reassign sequential IDs from 1 to table length"""
    cur.execute(f"""
        WITH numbered AS (
            SELECT ctid, ROW_NUMBER() OVER (ORDER BY id) as new_id FROM {table}
        )
        UPDATE {table} SET id = numbered.new_id FROM numbered WHERE {table}.ctid = numbered.ctid
    """)
    conn.commit()
    cur.execute(f"SELECT COUNT(*) FROM {table}")
    count = cur.fetchone()[0]
    print(f"[REINDEX] {table}: Reindexed ids from 1 to {count}")

def impute_tables(nStocksLog):
    """Main imputation procedure for a given nStocksLog value"""
    conn, cur = get_conn_cur(config)

    validate_table_name = f"stocks_{nStocksLog}_validate"
    prefix = f"stocks_{nStocksLog}_"

    print(f"\n=== Finding violating IDs from {validate_table_name} ===")
    violating_ids = get_violating_ids(cur, validate_table_name)
    print(f"Found {len(violating_ids)} violating IDs with std > 6500")

    if len(violating_ids) == 0:
        print("No violating IDs found. Nothing to impute.")
        cur.close()
        conn.close()
        return

    print(f"\n=== Getting all tables matching {prefix}* ===")
    related_tables = get_related_tables(cur, prefix)
    # Remove validate table from list - we'll process it last
    other_tables = [t for t in related_tables if t != validate_table_name]
    print(f"Found tables: {other_tables}")
    print(f"(Validate table {validate_table_name} will be processed last)")

    print(f"\n=== Filtering tables ===")
    for table in other_tables:
        filter_table(conn, cur, table, violating_ids)
    # Filter validate table last
    filter_table(conn, cur, validate_table_name, violating_ids)

    print(f"\n=== Reindexing tables ===")
    for table in other_tables:
        reindex_table(conn, cur, table)
    # Reindex validate table last
    reindex_table(conn, cur, validate_table_name)

    print(f"\n=== Imputation complete ===")
    cur.close()
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python impute.py <nStocksLog>")
        print("Example: python impute.py 5")
        sys.exit(1)

    nStocksLog = int(sys.argv[1])
    impute_tables(nStocksLog)
