#!/usr/bin/python3
import argparse
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from joblib import delayed,Parallel
import os
parser = argparse.ArgumentParser(description='Maximum likelihood inference of bursting kinetics from scRNA-seq data')
parser.add_argument('file', metavar='file', type=str, nargs=1,help='.pkl file with allelic-resolution transcript counts' )
args = parser.parse_args()
filename = args.file[0]

outfile = filename + '.csv'
print('Saving result to ' + outfile)
df = pd.read_pickle(filename)
df.to_csv(outfile)

