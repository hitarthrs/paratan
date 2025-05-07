#!/openmc_venv/bin/python
#
# for python3
#
import argparse
import pandas as pd
import numpy as np
#
parser = argparse.ArgumentParser(description=(
             'Reads a surface source data file and bins Zo using numpy.'
             'Also calculates average Ro in each Zo bin using pandas dataframe.'
             'Assumes surface source file format is Xo Yo Zo.'
             'Surface source file must have extra spaces removed.'))
# data file needs to be clean (remove excess spaces from file)
# cat surface_source.txt | tr -s ' ' | sed "s/^ //g" | sed "s/ $//g">surface_sourceClean.txt
#
parser.add_argument('filename1', help='Name of data file to read.')
parser.add_argument('-b', dest='numbins', help=('Number of bins.'),type=int)
#
parser.add_argument('rangemin', help='Range minimum.', type=float)
parser.add_argument('rangemax', help='Range maximum.', type=float)
#
args = parser.parse_args()
#
#
filein1=args.filename1
num_databins=args.numbins if args.numbins is not None else 10
range_min=args.rangemin
range_max=args.rangemax
#
print("\n The inputfile is:      ", filein1, "\n")
print("The number of bins is: ",num_databins," Range is:",range_min, " to ", range_max, "\n")
#
f1=open(filein1,'r')
#
# let numpy histogram function set-up binEdges array and bin the data
#
# read in Zo values from data file to bin with numpy (easier than with pandas)
data_array=np.genfromtxt(filein1,usecols=(2)) # read in 3rd column (Zo)
# generate the binned data using histogram function
(binCounts,binEdges)=np.histogram(data_array, bins=num_databins, range=(range_min,range_max),density=False)
print("Type for binCounts is: ", type(binCounts), " Shape is: ",binCounts.shape)
print("Type for binEdges is: ", type(binEdges), " Shape is: ",binEdges.shape)
#
print("Bin edges: \n",binEdges) # boundaries of each bin
print("Bin counts: \n",binCounts) # counts in each bin
#
np.savetxt('testbinEdges.txt', binEdges,fmt="%.3e")
np.savetxt('testbinCounts.txt', binCounts,fmt="%.3e")
#
########################################################################
# now calculate average Ro in each bin using pandas dataframe
print("\n Reading in datafile to pandas dataframe...")
# read in data file to a pandas dataframe with space as delimiter (need to determine average Ro for each bin/bucket
df=pd.read_csv(filein1, sep=' ', header=None, names=["Xo", "Yo", "Zo"])
# data file needs to be clean (remove excess spaces from file)
# cat surface_source.txt | tr -s ' ' | sed "s/^ //g" | sed "s/ $//g">surface_sourceClean.txt
#
print("The dataframe is: \n", df)
#
# calculate Ro and add to the dataframe
df['Ro']=((df['Xo']*df['Xo'])+(df['Yo']*df['Yo']))**0.5
print("The dataframe with Ro is: \n", df)
#
# bin or bucket the data
out,bucket_boundaries=pd.cut(df.Zo, bins=binEdges, include_lowest=True, right=False, retbins=True)
bucket_counts=out.value_counts()
df["bucket"] = pd.cut(df.Zo, binEdges) # add bucket column to df for later processing, hard code in Zo
print("The dataframe with buckets is: \n", df)
#
# check that data frame binning/bucketing matches numpy
print("\n Type for bucket_boundaries is: ", type(bucket_boundaries), " Shape is: ",bucket_boundaries.shape)
print("\n bin boundaries \n", bucket_boundaries)
print("\n Type for bucket_counts is: ", type(bucket_counts), " Shape is: ",bucket_counts.shape)
print("\n counts \n", bucket_counts) # note prints by number of counts in descending order
#
print("\n Next: Group by bin/bucket and aggregate Ro with mean/average: \n")
dfRo=df.groupby('bucket').agg(
    Ro_mean=pd.NamedAgg(column="Ro", aggfunc="mean")
    )
print(dfRo)
#
dfRo.to_csv('testromean.txt', sep=' ', header=False) # writes whole dataframe: bucket and ro_mean
#
#
f1.close()
#
