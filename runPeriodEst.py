import os
import sys
from fourierSeries import estSigGStat
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

p,w,g = estSigGStat(x,t,wV = wV, nPerm=perms,permEst=False,damp=damp)
if np.isnan(p) or p==0:
	p,w,g = estSigGStat(x,t,wV = wV, nPerm=perms*100,permEst=False,damp=damp)

#f = open(outPath,'a')
#f.write('>>\t'+str(index)+'\t'+str(p)+'\t'+str(w)+'\n')

#f.flush()
#os.fsync(f)
#f.close()
if ~np.isnan(p): out = '>>\t'+str(index)+'\t%05.4E\t%05.4E\t%05.4E\n' % (p,w,g)
else: out = '>>\t'+str(index)+'\tnan\tnan\tnan\n'

print(out)

