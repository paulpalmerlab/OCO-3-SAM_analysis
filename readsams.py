import numpy as np
import pandas as pd
import sys
import os

filename = 'SAM_TG.csv'

p = pd.read_csv(filename,header=0)

p = np.array(p)

#getfiles = []

nsams = len(p[:,0])

for ii in np.arange(nsams):
    getfilename = p[ii,4]+p[ii,3]
    print(getfilename)
    #getfiles.append(p[ii,4]+p[ii,3])
    os.system('scp ppalmer@oil.jpl.nasa.gov:'+getfilename+' .')
    
