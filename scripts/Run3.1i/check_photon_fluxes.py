import os
import glob
import sqlite3
import subprocess
import numpy as np
import pandas as pd

#db_tables = {
#    'agn_truth_cat.db': 'agn_variability_truth',
#    'lensed_agn_truth_cat.db': 'lensed_agn_variability_truth',
#    'lensed_sne_truth_cat.db': 'lensed_sn_variability_truth'
#}

prod_dir = '/global/cscratch1/sd/descim/Run3.1i/truth_cats'
db_tables = {
    os.path.join(prod_dir, 'agn_variability_truth_cat.db'):
    'agn_variability_truth',
    os.path.join(prod_dir, 'lensed_agn_variability_truth_cat.db'):
    'lensed_agn_variability_truth',
    os.path.join(prod_dir, 'lensed_sne_truth_cat.db'):
    'lensed_sn_variability_truth'
}

sims_dir = '/global/cfs/cdirs/lsst/production/DC2_ImSim/Run3.1i/sim/y?-ddf'

for db_file, table in db_tables.items():
    print(db_file)
    with sqlite3.connect(db_file) as con:
        df = pd.read_sql(f'''select id, obsHistID, num_photons from {table}
                             order by RANDOM() limit 100''', con)

    for _, row in df.iterrows():
        visit = row['obsHistID']
        object_id = row['id']
        try:
            pattern = os.path.join(sims_dir, f'{visit:08d}')
            visit_dir = glob.glob(pattern)[0]
        except IndexError:
            continue
        command = f'gunzip -c {visit_dir}/centroid*.gz | grep {object_id}'
        try:
            results = subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError:
            continue
        photon_flux = float(results.split()[1])
        dflux = row['num_photons'] - photon_flux
        print(object_id, row['num_photons'], photon_flux, dflux)
        if np.abs(dflux) > 1:
            print(results)
    print()

db_file = os.path.join(prod_dir, 'lensed_host_truth_cat.db')
table = 'truth_summary'
print(db_file)
with sqlite3.connect(db_file) as con:
    df = pd.read_sql(f'select * from {table}', con)

visit_dirs = glob.glob(os.path.join(sims_dir, '0*'))
np.random.shuffle(visit_dirs)

for visit_dir in visit_dirs[:6]:
    band = glob.glob(os.path.join(visit_dir, 'lsst*.fits'))[0]\
           [-len('x.fits'):-len('.fits')]
    for iloc in np.random.choice(range(len(df)), 6):
        row = df.iloc[iloc]
        object_id = row['id']
        num_photons = row[f'num_photons_{band}']
        command = f'gunzip -c {visit_dir}/centroid*.gz | grep {object_id}'
        try:
            results = subprocess.check_output(command, shell=True)\
                                .decode('utf-8').strip().split('\n')
        except subprocess.CalledProcessError:
            continue
        xpix, ypix = [float(_) for _ in results[0].split()[3:5]]
        if xpix < 0 or xpix > 4000 or ypix < 0 or ypix > 4000:
            #print("skipping", results)
            continue
        photon_flux = sum([float(_.split()[1]) for _ in results])
        dflux = num_photons - photon_flux
        print(object_id, num_photons, photon_flux, dflux)
        if np.abs(dflux) > 1:
            print('\n'.join(results))
