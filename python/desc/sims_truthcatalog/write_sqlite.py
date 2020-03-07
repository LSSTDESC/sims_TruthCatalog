"""
Function to write truth_summary table
"""
import sqlite3
from . import sqlite_utils


__all__ = ['write_sqlite']


def write_sqlite(dbfile, ids, galaxy_ids, ra, dec, redshift, is_variable,
                 is_pointsource, flux_by_band_MW, flux_by_band_noMW, good_ixes):
    """
    Write the truth_summary table.
    """
    with sqlite3.connect(dbfile) as conn:
        sqlite_utils.write_column_descriptions(conn)
        cursor = conn.cursor()

        cmd = '''CREATE TABLE IF NOT EXISTS truth_summary
        (id TEXT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
        redshift FLOAT, is_variable INT, is_pointsource INT,
        flux_u FLOAT, flux_g FLOAT, flux_r FLOAT,
        flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
        flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT,
        flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT)'''
        cursor.execute(cmd)
        conn.commit()

        values = ((str(ids[i_obj]), int(galaxy_ids[i_obj]),
                   ra[i_obj], dec[i_obj], redshift[i_obj],
                   is_variable[i_obj], is_pointsource[i_obj],
                   flux_by_band_MW['u'][i_obj], flux_by_band_MW['g'][i_obj],
                   flux_by_band_MW['r'][i_obj], flux_by_band_MW['i'][i_obj],
                   flux_by_band_MW['z'][i_obj], flux_by_band_MW['y'][i_obj],
                   flux_by_band_noMW['u'][i_obj], flux_by_band_noMW['g'][i_obj],
                   flux_by_band_noMW['r'][i_obj], flux_by_band_noMW['i'][i_obj],
                   flux_by_band_noMW['z'][i_obj], flux_by_band_noMW['y'][i_obj])
                  for i_obj in good_ixes)

        cursor.executemany('''INSERT INTO truth_summary
                              VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                           values)
        conn.commit()

        # index to speed up location searches
        cursor.execute('create index radec_ix on truth_summary(ra,dec)')
        conn.commit()
