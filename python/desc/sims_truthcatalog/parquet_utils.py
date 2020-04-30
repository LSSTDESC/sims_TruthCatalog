import os
import sys

import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd
import sqlite3

from desc.sims_truthcatalog.script_utils import print_callinfo

__all__ = ["convert_sqlite_to_parquet", "compare_sqlite_parquet"]
def  _transpose(records, column_dict, schema, n_rec=None):
    '''
    return data as list of columns
    '''
    data_dict = {k : [] for k in column_dict.keys()}
    dat = list(data_dict.values())
    rng = len(records)
    if n_rec is not None: rng = n_rec
    ic = len(column_dict.keys())
    for ir in range(rng):
        for i in range(ic):
            dat[i].append(records[ir][i])

    print("Type of data ", type(dat))
    print("len of data ", len(dat))
    for c in dat:
        print("An item: ", c[0], " its type: ", type(c[0]))
    print(data_dict['id'][0])

    return pa.Table.from_arrays(dat, schema=schema)

def convert_sqlite_to_parquet(dbfile, pqfile, table,
                              n_group=None, max_group_gbyte=1.0, order_by = None):
    '''
    Write a parquet file corresponding to contents of a table from an sqlite3 db.

    Parameters:
    dbfile          input sqlite3 
    pqfile          path for output file
    table           write contents of this table to parquet
    n_group         # of row groups to write in output parquet file. If not 
                    specified, use max_group_gbyte to determine workable value
    max_group_gbyte maximum size for a row group.  Will write more row groups
                    than specified by ngroup if need be.
    '''

    statinfo = os.stat(dbfile)
    min_groups = np.ceil(statinfo.st_size / (float(max_group_gbyte) * 10e9))
    ng = min_groups
    if n_group is not None: ng = max(min_groups, n_group)

    column_dict = {}
    with sqlite3.connect(dbfile) as conn:
        cursor = conn.cursor()
        cmd = 'select count(*) from ' + table
        res = cursor.execute(cmd)
        total_row = res.fetchone()[0]
        row_per_group = total_row
        if ng > 1:
            row_per_group = int(np.ceil(total_row/float(ng)))

        # Get names, types
        meta_res = cursor.execute('PRAGMA table_info({})'.format(table))
        column_dict = {t[1]: t[2] for t in meta_res.fetchall()}

        type_translate = {'BIGINT' : 'int64', 'INT' : 'int32', 'INTEGER' : 'int32',
                          'FLOAT' : 'float32', 'DOUBLE' : 'float64', 'TEXT' : 'string'}
        type_pa = {'int64' : pa.int64(), 'int32' : pa.int32(),
                   'float32': pa.float32(), 'float64' : pa.float64(), 'string' : pa.string()}

        for (k,v) in column_dict.items():
            if v in type_translate.keys():
                column_dict[k] = type_translate[v]
            else:
                print(f"For key {k} found unknown type {v}, setting to float32")
                column_dict[k] = 'float32'

        print('column_dict: ')
        for (k, v) in column_dict.items():
            print(k, " : ", v, " : ")
        # Make the parquet schema
        fields = []
        for (k, v) in column_dict.items():
            fields.append((k, type_pa[v]))
        schema = pa.schema(fields)
        for k in schema:
            print(k)
                
        done = False

        # Fetch the sqlite data

        prev_row = 0
        limit_row = min(row_per_group, total_row)

        writer = pq.ParquetWriter(pqfile, schema)

        while limit_row <= total_row:
            cmd = f"select * from {table} where rowid > {prev_row} and rowid <= {limit_row}"
            
            if order_by is not None:
                cmd += ' order by ' + order_by

            print('\nFetch command is:\n', cmd)
            
            cursor.execute(cmd)

            #with pq.ParquetWriter(pqfile, schema) as writer:
            #    #while not done:
            #    #records =  cursor.fetchmany(row_per_group)
            records =  cursor.fetchall()
            if len(records) == 0:
                done = True
                break
            to_write = _transpose(records, column_dict, schema)
            writer.write_table(to_write)
            if len(records) < row_per_group:
                done = True
                break
            prev_row = limit_row
            limit_row = min(prev_row + row_per_group, total_row)
    print('Conversion successful')                   #  ***DEBUG***
        
def compare_sqlite_parquet(sqlite_file, parquet_file, sqlite_table, id_column=None, n_rows=100,
                           verbose=False, check_cols=None):
    '''
    Compare an sqlite table with one from a parquet file

    Parameters
    ----------
    sqlite_file     string    path to sqlite file
    parquet_file    string    path to parquet file
    sqlite_table    string    sqlite table name
    id_column       string    If not none, a unique id which exists in both tables
    n_rows          int       If there is an id_column, number of rows to compare directly;
                              else ignored
    check_cols      list      compare values for these columns

    Return
    ------
    True if checks pass
    False if there is a failure.  In this case raise warning with message
    '''

    # Make sure we can connect to both files
    try:
        sq_conn = sqlite3.connect(sqlite_file)
    except Exception as ex:
        warnings.warn('Failed to connect to sqlite file with error {}'.format(str(ex)))
        return False

    sq_conn.row_factory = sqlite3.Row
    sq_cur = sq_conn.cursor()
    sqlite_desc_query = 'select * from {} limit 1'.format(sqlite_table)
    sq_cur.execute(sqlite_desc_query)

    if verbose: print('Opened sqlite; executed one-row query')

    try:
        pq_file = pq.ParquetFile(parquet_file)
    except Exception as ex:
        warnings.warn('Failed to open parquet file with error {}'.format(str(ex)))
        return False

    if verbose: print('Opened parquet file')
    
    # Check that schemas are compatible (i.e., column names match.  Ideally should
    # also check types are compatible)
    pq_lst = []
    if verbose: print(pq_file.schema)
    for i in range(len(pq_file.schema)):
        nm = str(pq_file.schema[i].name)
        pq_lst.append(nm)
    pq_set = set(pq_lst)

        
    #sql_set_lc = set()
    sql_lst = []
    for n in sq_cur.description:
        sql_lst.append(str(n[0]))
    sql_set = set(sql_lst)
    if verbose:
        print('sql list has {} members'.format(len(sql_lst)))
        for m in sql_lst:
            print(m)
        print('pq list has {} members'.format(len(pq_lst)))
        for m in pq_lst:
            print(m)

    if sql_set != pq_set:
        warnings.warn('Table columns do not match')
        return False

    # Check that number of rows match
    cnt_q = 'select count({}) from {}'.format(sql_lst[0],sqlite_table)
    sq_cur.execute(cnt_q)
    sq_num = sq_cur.fetchone()[0]

    if sq_num != pq_file.metadata.num_rows:
        warnings.warn('sqlite row cnt {} != parquet row cnt {}'.format(sq_num,
                                                                       pq_file.metadata.num_rows))
        return False

    # If there is an id_colum    read in first n_rows from sqlite file
    #     and for each row  attempt to select the 'same' rows from parquet; compare

    if id_column is not None:
        sq_query = 'select * from {sqlite_table} limit {n_rows}'.format(**locals())
        sq_cur.execute(sq_query)
        df = pd.read_parquet(parquet_file, engine='pyarrow')
        
        r = sq_cur.fetchone()
        while r is not None:
            pq_row = df.loc[lambda d: d[id_column] == r[id_column], :]
            if pq_row.shape[0] == 0:
                warnings.warn('sqlite row missing in parquet')
                return False
            if pq_row.shape[0] > 1 :
                warnings.warn('id column values not unique in parquet')
                return False
            r = sq_cur.fetchone()
            

    return True
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Write parquet with content from sqlite or compare existing files')
    parser.add_argument('dbfile', type=str, help='input sqlite file; required')
    parser.add_argument('--pqfile', type=str,
                        help='parquet filepath. Defaults to test.parquet if output')
    parser.add_argument('--table', type=str, help='table to be written or compared parquet',
                        default='truth_summary')
    parser.add_argument('--n-group', type=int, default=None,
                        help='number of row groups. By default compute from max_group_gbyte')
    parser.add_argument('--max-group-gbyte', type=float, default=1.0,
                        help='max allowed size in gbytes of a row group')
    parser.add_argument('--check', action='store_true',
                        help='If set, compare sqlite and parquet')
    parser.add_argument('--id-column', default=None)
    parser.add_argument('--n-check', type=int, default=10,
                        help='Number of rows from sqlite to check. Ignored in no check or id_colume is None')

    args = parser.parse_args()

    print_callinfo(sys.argv[0], args)
    #announce='\nCalled with arguments\ndbfile={}\npqfile={}\ntable={} n_groups={}\n'
    #announce += 'max_group_gbyte={} '
    #print(announce.format(args.dbfile, args.pqfile, args.table, args.n_group,
    #                      args.max_group_gbyte))

    if not args.check:
        convert_sqlite_to_parquet(args.dbfile, args.pqfile, args.table,
                                  n_group=args.n_group,
                                  max_group_gbyte=args.max_group_gbyte, order_by = 'rowid')
    else:
        ok = compare_sqlite_parquet(args.dbfile, args.pqfile, args.table, args.id_column,
                                    args.n_check)
        if ok: print('Comparison succeeded')

