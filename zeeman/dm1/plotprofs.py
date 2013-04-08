#!/usr/bin/python
import numpy as np
import pylab as pl
import matplotlib.pyplot as mpt
import sys

b=sys.argv[1]
px,val = pl.loadtxt(b+'.txt',unpack=True,usecols=[0,1])
mpt.figure(figsize=(13,0.7),dpi=150)
mpt.plot(px,val-20,'b,')
mpt.title('b'+str(b))
mpt.subplots_adjust(left=0,right=1,top=0.9,bottom=0.1)
if b=='225':
    mpt.savefig('plot225.png')
mpt.show()
raw_input()
