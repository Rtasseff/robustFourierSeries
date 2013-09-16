import os
import sys
from fourierSeries import estSigFreqInd
import numpy as np

index = int(sys.argv[1])
tPath = sys.argv[2]
xPath = sys.argv[3]
wPath = sys.argv[4]
perms = int(sys.argv[5])
damp = float(sys.argv[6])
#outPath = sys.argv[7]

	
t = np.loadtxt(tPath)
X = np.loadtxt(xPath)
if wPath == 'NA':
	wV = []
else:
	wV = np.loadtxt(wPath)

x = X[index,:]

p,w = estSigFreqInd(x,t,wV, nPerm=perms,damp=damp)
if np.any(np.isnan(p)) or np.any(p==0):
	p,w = estSigFreqInd(x,t,wV, nPerm=perms*100,damp=damp)

#f = open(outPath,'a')
#f.write('>>\t'+str(index)+'\t'+str(p)+'\t'+str(w)+'\n')

#f.flush()
#os.fsync(f)
#f.close()

out = '>>\t'+str(index)

for i in range(len(p)):
	if ~np.isnan(p[i]): out = '%s\t%05.4E' % (out,p[i])
	else: out = out+'\tnan'

out = out+'\n'

print out
