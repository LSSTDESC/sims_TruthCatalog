
from datetime import datetime as dt
import shutil
import os
import sqlite3

__all__ = ["write_column_descriptions", "merge_sqlite3_dbs"]

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


def merge_sqlite3_dbs(infiles, outfile):
    """
    Function to merge the db tables in a list of input sqlite3 files.
    The corresponding tables in each file are assumed to have the same
    schemas.

    Parameters
    ----------
    infiles: list-like
        List of sqlite3 filenames.
    outfile: str
        Filename of output sqlite3 file to contain the merged tables.
    """
    shutil.copy(infiles[0], outfile)
    with sqlite3.connect(outfile) as conn:
        for infile in infiles[1:]:
            conn.execute(f"attach '{infile}' as in_db")
            for row in conn.execute("""select * from in_db.sqlite_master
                                       where type='table'"""):
                combine = f"insert into {row[1]} select * from in_db.{row[1]}"
                conn.execute(combine)
            conn.commit()
            conn.execute("detach database in_db")
