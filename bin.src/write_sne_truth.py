#!/usr/bin/env python
"""
Script to write the SNe truth catalog.
"""
import argparse
import desc.sims_truthcatalog

parser = argparse.ArgumentParser(description="Write a truth catalog (sqlite) for SNe given a SNe parameters db file")
parser.add_argument('outfile', type=str,
                    help='Filename for output sqlite file')
parser.add_argument('sne_db_file', type=str,
                    help='Filename for input SNe parameters db file')

args = parser.parse_args()

writer = desc.sims_truthcatalog.SNeTruthWriter(args.outfile, args.sne_db_file)
writer.write()
