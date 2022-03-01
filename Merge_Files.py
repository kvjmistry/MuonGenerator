import pandas as pd
import sys
import glob

filewildcard = sys.argv[1]
files = glob.glob(filewildcard)

print(files)

Data = []

for f in files:
    Data.append(pd.read_hdf(f,'Data'))

Data_m = pd.concat(Data)

print(Data_m)

file_out = sys.argv[2]
Data_m.to_hdf(file_out,'Data')