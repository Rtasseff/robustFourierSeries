import numpy as np

	
listId = 'periodogramOut_20130301_Lin2009_log2_periodic_id.dat'
fin = open (listId,'r')
x = np.loadtxt('X.dat')
n,m = x.shape
w = np.loadtxt('w.dat')
m = len(w)
p = np.ones(n)*2
w = np.zeros(n)
missedOut = open('missed_in.txt','w')

for line in fin:
	data = line.split()
	if len(data)>0 and data[0]=='>>':
		index = int(data[1])
		p[index] = float(data[2])
		w[index] = float(data[3])

missed = p>1
for i in range(len(missed)):
	if missed[i]:
		missedOut.write(str(i)+'\n')

missedOut.close()
np.savetxt('p.dat',p)
np.savetxt('wMax.dat',w)
	
	

