#!/usr/bin/env python
# coding: utf-8

import configparser
import psycopg2
from itertools import combinations

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

def get_tables_by_prefix(cur, prefix):
    """Get all tables matching a prefix"""
    cur.execute("""
        SELECT table_name FROM information_schema.tables
        WHERE table_schema = 'public' AND table_name LIKE %s
        ORDER BY table_name
    """, (f"{prefix}%",))
    return [row[0] for row in cur.fetchall()]

def compare_two_tables(cur, table1, table2):
    """Compare id and price columns between two tables"""
    # Check if they have the same number of rows
    cur.execute(f"SELECT COUNT(*) FROM {table1}")
    count1 = cur.fetchone()[0]
    cur.execute(f"SELECT COUNT(*) FROM {table2}")
    count2 = cur.fetchone()[0]

    if count1 != count2:
        return False, f"Row count mismatch: {table1}={count1}, {table2}={count2}"

    # Check if all (id, price) pairs match
    cur.execute(f"""
        SELECT COUNT(*) FROM (
            SELECT id, price FROM {table1}
            EXCEPT
            SELECT id, price FROM {table2}
        ) diff
    """)
    diff_count = cur.fetchone()[0]

    if diff_count > 0:
        return False, f"{diff_count} rows in {table1} not found in {table2}"

    # Check reverse direction
    cur.execute(f"""
        SELECT COUNT(*) FROM (
            SELECT id, price FROM {table2}
            EXCEPT
            SELECT id, price FROM {table1}
        ) diff
    """)
    diff_count_rev = cur.fetchone()[0]

    if diff_count_rev > 0:
        return False, f"{diff_count_rev} rows in {table2} not found in {table1}"

    return True, "OK"

def verify_all_pairs(prefix):
    """Verify all pairs of tables matching the prefix have identical id and price"""
    conn, cur = get_conn_cur(config)

    tables = get_tables_by_prefix(cur, prefix)
    print(f"Found {len(tables)} tables matching '{prefix}*':")
    for t in tables:
        print(f"  - {t}")

    if len(tables) < 2:
        print("Need at least 2 tables to compare")
        return

    print(f"\n=== Comparing all pairs ===")

    # Use first table as reference instead of all pairs (more efficient)
    reference_table = tables[0]
    all_match = True

    print(f"\nUsing {reference_table} as reference table:")
    for table in tables[1:]:
        match, msg = compare_two_tables(cur, reference_table, table)
        status = "✓" if match else "✗"
        print(f"  {status} {reference_table} vs {table}: {msg}")
        if not match:
            all_match = False

    if all_match:
        print(f"\n✓ All tables have identical (id, price) values!")
    else:
        print(f"\n✗ Some tables have mismatched (id, price) values!")

    cur.close()
    conn.close()

def verify_specific_pair(table1, table2):
    """Verify a specific pair of tables"""
    conn, cur = get_conn_cur(config)

    match, msg = compare_two_tables(cur, table1, table2)
    status = "✓" if match else "✗"
    print(f"{status} {table1} vs {table2}: {msg}")

    if not match:
        # Show some example differences
        print("\nExample differences:")
        cur.execute(f"""
            (SELECT 'in {table1} only' as source, id, price FROM {table1}
             EXCEPT
             SELECT 'in {table1} only', id, price FROM {table2})
            UNION ALL
            (SELECT 'in {table2} only' as source, id, price FROM {table2}
             EXCEPT
             SELECT 'in {table2} only', id, price FROM {table1})
            LIMIT 10
        """)
        for row in cur.fetchall():
            print(f"  {row}")

    cur.close()
    conn.close()

if __name__ == "__main__":
    import sys

    if len(sys.argv) == 2:
        # Verify all tables with a prefix
        prefix = sys.argv[1]
        verify_all_pairs(prefix)
    elif len(sys.argv) == 3:
        # Verify specific pair
        verify_specific_pair(sys.argv[1], sys.argv[2])
    else:
        print("Usage:")
        print("  python verify_tables.py <prefix>              # Compare all tables with prefix")
        print("  python verify_tables.py <table1> <table2>     # Compare two specific tables")
        print("")
        print("Examples:")
        print("  python verify_tables.py stocks_4_")
        print("  python verify_tables.py stocks_4_10_seeded_1 stocks_4_10_seeded_2")