import numpy as np
import fourierSeries


t = np.array([8.67, 10.52, 12.42, 15.67, 20.27, 23.52, 27.52, 30.52, 33.77, 37.32, 40.67, 45.22, 48.52, 51.67, 55.17, 58.47, 62.52, 65.52])

t = np.arange(18)

n = len(t)

e = .25

x = np.zeros((200,n))

wV = (2*np.pi/n)*np.arange(1,int(n/2))

for i in range(100):
	x[i,:] =  e*np.random.randn(n)

b = np.sqrt(2)
for i in range(100,200):
	w = np.random.rand(1)*(np.max(wV)-np.min(wV))+np.min(wV)
	s = np.random.rand(1)*2.0*n - n
	x[i,:] = b*np.cos(w*t+s) + e*np.random.randn(n)

#g = np.zeros(200)
#for i in range(200):
#	g[i],wTmp = fourierSeries.gStatMax(x[i,:],t,wV)
#
#np.savetxt('g.dat',g)

p = np.zeros(200)
for i in range(200):
	p[i],wTmp = fourierSeries.estSigGStat(x[i,:],t,wV,300,False)

np.savetxt('p.dat',p)

tmp = -1*np.log10(p)
m = np.median(tmp)
plt.plot(tmp,'o');plt.plot([0,len(p)],[m,m]);plt.show()




