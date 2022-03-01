import pandas as pd
import sys
import glob
from invisible_cities.io.mcinfo_io import load_mcparticles_df

filewildcard = sys.argv[1]
files = glob.glob(filewildcard)

print(filewildcard)


data = load_mcparticles_df(filewildcard) 

print(data)

file_out = sys.argv[2]
pd.DataFrame(data).to_hdf(file_out,'Data')
