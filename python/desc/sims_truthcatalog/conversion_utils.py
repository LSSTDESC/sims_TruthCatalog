import sqlite3
import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd

from datetime import datetime as dt
import os

__all__ = ["convert_sqlite_to_parquet", "write_column_descriptions"]

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

        type_translate = {'BIGINT' : 'int64', 'INT' : 'int32', 'FLOAT' : 'float32',
                          'DOUBLE' : 'float64'}
        type_pa = {'int64' : pa.int64(), 'int32' : pa.int32(),
                   'float32': pa.float32(), 'float64' : pa.float64()}

        for (k,v) in column_dict.items():
            if v in type_translate.keys():
                column_dict[k] = type_translate[v]

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
        cmd = 'select * from ' + table
        if order_by is not None:
            cmd += ' order by ' + order_by

        print('\nFetch command is:\n', cmd)
            
        cursor.execute(cmd)

        with pq.ParquetWriter(pqfile, schema) as writer:
            while not done:
                records =  cursor.fetchmany(row_per_group)
                if len(records) == 0:
                    done = True
                    break
                to_write = _transpose(records, column_dict, schema)
                writer.write_table(to_write)
                if len(records) < row_per_group:
                    done = true
    print('Conversion successful')                   #  ***DEBUG***
        
def write_column_descriptions(conn):
    cursor = conn.cursor()
    try:
        cursor.execute('select count(*) from column_descriptions')
        return
    except:
        pass

    cmd = '''CREATE TABLE column_descriptions
    (name text, description text, dtype text)'''
    cursor.execute(cmd)
    conn.commit()
    rowdat = (('id', 'Unique identifier for object (star or galaxy)', 'int8'),
              ('host_galaxy', '-1 if no host', 'int8'),
              ('ra', 'Right ascension, lensed (deg)', 'float8'),
              ('dec', 'Declination, lensed (deg)', 'float8'),
              ('redshift', 'Same as CosmoDC2 redshift (not redshift_true)',
               'float4'),
              ('is_variable', '1 if variable; else 0', 'int4'),
              ('is_pointsource', '1 if star; 0 for galaxy', 'int4'),
              ('flux_u', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_g', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_r', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_i', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_z', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_y', 'Time-avged flux, lensed, all extinction incl, (nanoJ)',
               'float4'),
              ('flux_u_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
              ('flux_g_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
              ('flux_r_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
              ('flux_i_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
              ('flux_z_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
              ('flux_y_noMW', 'Time-avged flux, lensed, no MW extinct, (nanoJ)',
               'float4'),
    )

    for i in range(len(rowdat)):
        cursor.execute('''INSERT INTO column_descriptions VALUES (?,?,?)''',
                       (rowdat[i][0], rowdat[i][1], rowdat[i][2]))
        conn.commit()
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Write parquet with content from sqlite')
    parser.add_argument('dbfile', type=str, help='input sqlite file; required')
    parser.add_argument('--pqfile', type=str, help='output filepath',
                        default='test.parquet')
    parser.add_argument('--table', type=str, help='table to be written to parquet',
                        default='truth_summary')
    parser.add_argument('--n-group', type=int, default=None,
                        help='number of row groups. By default compute from max_group_gbyte')
    parser.add_argument('--max-group-gbyte', type=float, default=1.0,
                        help='max allowed size in gbytes of a row group')

    args = parser.parse_args()
    
    announce='\nCalled with arguments\ndbfile={}\npqfile={}\ntable={} n_groups={}\n'
    announce += 'max_group_gbyte={} '
    print(announce.format(args.dbfile, args.pqfile, args.table, args.n_group,
                          args.max_group_gbyte))

    convert_sqlite_to_parquet(args.dbfile, args.pqfile, args.table,
                              n_group=args.n_group,
                              max_group_gbyte=args.max_group_gbyte, order_by = 'id')
