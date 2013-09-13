import numpy as np
#listId = '20120413'
#name = 'WT'
#mypath='/titan/cancerregulome9/workspaces/pngisb/rtasseff/'
#f=open('golemJobList_'+listId+'_'+name+'.txt','w')
#for i in range(10715):
#	f.write('1 /hpc/bin/python2.7 '+mypath+'runPeriodEst.py '+ \
#		str(i)+' '+mypath+'t'+name+'.dat '+mypath+'x'+name+'.dat '+mypath+'w.dat 1000 '+mypath+'periodogramOut_'+listId+'_'+name+'.dat\n')
#
#f.close()

dampeningCeof = 0
label = '20130913'
name = 'GSS0373'
n =  45141
mypath='/titan/cancerregulome9/workspaces/pngisb/rtasseff/periodID5/'
f=open('golemJobList_'+label+'_'+name+'.txt','w')
#for i in range(n):
for i in range(n):
	f.write('1 /tools/bin/python2.7 '+mypath+'runPeriodEst.py '+ \
		str(i)+' '+mypath+'t.dat '+mypath+'X.dat '+mypath+'w.dat 1000 '+str(dampeningCeof)+' '+mypath+'periodogramOut_'+label+'_'+name+'.dat\n')

f.close()


