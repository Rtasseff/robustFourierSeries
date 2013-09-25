#
# 
#     Copyright (C) 2003-2012 Institute for Systems Biology
#                             Seattle, Washington, USA.
# 
#     This library is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 2.1 of the License, or (at your option) any later version.
# 
#     This library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with this library; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
#
#
######################
""" This module contains methods to calculate
the periodogram for time series data by
estimating the coefficents of a discreat 
Fourier Series.  The g-statistic is also
estimated and tested to identify the periodic 
component that is statistically significant.
This uses the robust non-uniform sampling 
algorithm described in Ahdesmai 2007
(Robust regression for periodicity detection in non-uniformly 
sampled time-course gene expression data) 
and the genral pereto distribution estimates
from Knijnenburg 2009
(Fewer permutations, more accurate P-values)

The basic method is estSigGStat; however,
many other methods are avalible.

Please excuse the spelling and poor commenting :)
"""


import numpy as np
# currently used for the robust regression:
import robust.robust_linear_model as sm
import robust.norms as norms
# currently used to improve pvalue estimates
import gpdPerm
# Created 201205 RAT
# Edited 20121205 RAT - added dampining term 





def estSigGStat(yPrime,tPrime,wV=[],nPerm=1000, permEst=False,damp=0,transform=False):
	"""Estimates the statistical significance 
	that the time series data has a periodic component
	corresponding to at least one of the frequencies 
	considered.  The periodogram is estimated via
	robust linear regression, the g statistic is used as the test
	statistic, permutation testing is done to estimate
	the null distribution (observations not related to time),
	and the general pereto distribution is used to improve
	estimation of the tail.
	Requires gpdPerm module
	y	array of data
	t	array of times corresponding to y
	nPerm	number of permutations to estimate null
	wV	array user defined frequencies
		if empty the time range and number of points
		will be used to estimate possible frequencies
	damp 	float, a priori damppining coef to be applied to data, exp{-damp*t}
	returns:
	p	scalar float- p value, probability of null
	w	scalar float- most weighted frequency (~largest periodic component) 
	permonly	if true will return standard permutation p statistic
			when gdpPerm fails, other wise nan will be returned
	"""
	import gpdPerm
	

	# check for nan's (its happend before)
	y,t = removenan(yPrime,tPrime)
	
	# get max
	g,w = gStatMax(y,t,wV,damp)
	gPerm = np.zeros(nPerm)
	# do the permutations
	for i in range(nPerm):
		yPerm = np.random.permutation(y)
		gPerm[i], _ = gStatMax(yPerm,t,wV,damp)
	if transform:p = gpdPerm.est(_transform(g),_transform(gPerm))
	else:p = gpdPerm.est(g,gPerm)
	if p==0:p = np.nan 
	if np.isnan(p) and permEst:
		p = (np.sum(gPerm>=g)+5.)/len(gPerm)

	return p,w,g
		
def _transform(t):
	return(1/(1-t))

def estSigFreqInd(y,t,w=[],nPerm=1000,damp=0):
	"""Estimates the statistical significance 
	of each specified frequancy.  The periodogram is estimated via
	robust linear regression, sum of squared trig coeff are used as the test
	statistic, permutation testing is done to estimate
	the null distribution (observations not related to time),
	and the general pereto distribution is used to improve
	estimation of the tail.
	Requires gpdPerm module
	y	array of data
	t	array of times corresponding to y
	nPerm	number of permutations to estimate null
	w	array user defined frequencies
		if empty the time range and number of points
		will be used to estimate possible frequencies
	damp	float, used as an a priori dampining coef e^{-damp*t} 
	returns:
	p	scalar float- p value, probability of null
	w	scalar float- most weighted frequency (~largest periodic component) 
	permonly	if true will return standard permutation p statistic
			when gdpPerm fails, other wise nan will be returned
	"""
	w = np.array(w)
	if len(w)==0:
		w = defultFreq(t)
	
	# check for nan's (its happend before)
	y,t = removenan(y,t)
	
	# get max
	B,b = calcPeriodogram(y,t,w,False,damp)	
	BPerm = np.zeros((len(B),nPerm))
	# do the permutations
	for i in range(nPerm):
		yPerm = np.random.permutation(y)
		BPermTemp, _ = calcPeriodogram(yPerm,t,w,False,damp)
		BPerm[:,i] = BPermTemp		
	
	p = _calcP(B,BPerm,nPerm)

	return p,w

def estSigFreq(y,t,w=[],nPerm=1000,permEst=False,damp=0):
	"""Estimates the statistical significance 
	of each specified frequancy.  The periodogram is estimated via
	robust linear regression, sum of squared trig coeff are used as the test
	statistic, permutation testing is done to estimate
	the null distribution (observations not related to time),
	and the general pereto distribution is used to improve
	estimation of the tail.
	Requires gpdPerm module
	y	array of data
	t	array of times corresponding to y
	nPerm	number of permutations to estimate null
	w	array user defined frequencies
		if empty the time range and number of points
		will be used to estimate possible frequencies
	damp	float, a priori dampining coef to be applied to data, exp^{-damp*t}
	returns:
	p	scalar float- p value, probability of null
	w	scalar float- most weighted frequency (~largest periodic component) 
	permonly	if true will return standard permutation p statistic
			when gdpPerm fails, other wise nan will be returned
	"""
	w = np.array(w)
	if len(w)==0:
		w = defultFreq(t)
	
	# check for nan's (its happend before)
	y,t = removenan(y,t)
	
	# get max
	B,b = calcPeriodogram2(y,t,w,damp)	
	BPerm = np.zeros((len(B),nPerm))
	# do the permutations
	for i in range(nPerm):
		yPerm = np.random.permutation(y)
		BPermTemp, _ = calcPeriodogram2(yPerm,t,w,damp)
		BPerm[:,i] = BPermTemp		
	
	p = _calcP(B,BPerm,nPerm,permEst)
	return p,w

def _calcP(stat,statPerm,nPerm, permEst=True):
	p = np.ones(len(stat))
	for i in range(len(stat)):
		p[i] = gpdPerm.est(stat[i],statPerm[i,:])
		if p[i]==0:p[i] = np.nan
		if np.isnan(p[i]) and permEst:
			p[i] = (np.sum(statPerm[i,:]>=stat[i])+5.)/nPerm
	return(p)



def gStatMax(y,t,wV=[],damp=0):
	"""calculates the maximum cycle g statistic for 
	the time series data provided. The statistic 
	describes the extent to which this time series 
	has at least one periodic component.
	Fits each possible frequency individually 
	and removes the residual for the next fit.
	The fit is performed twice; first the frequencies 
	are considered in the initial order, these estimates 
	are used to rank order the second round of fitting. 
	y	array - data vector
	t	array - time points corresponding to y
	damp	a priori dampping coef
	returns:
	g	float - the maximum g test statistic
	w	float - the frequency corresponding to g 
	"""
	if len(wV) == 0:
		wV = defultFreq(t)


	# calculate the periodogram 
	B,b = calcPeriodogram2(y,t,wV,damp)
	# find the maximum
	val = np.max(B)
	ind = np.argmax(B)
	# g stat
	
	g = val/np.sum(B)
	w = wV[ind]
	return g,w


def defultFreq(t):
	"""Since we are not using uniform time points
	the standard harmonic frequancies do not make sense
	we have, somewhat arbitrarily, chosen the following
	w_min corrisponds to a maximum period 
	which contains all points will be the range of 
	max(t) - min(t), and the w_max corrisponds to a 
	minimum period that contains at worst 2 points, 
	on a sorrted t -> t_sort, max_n(t_n+1-t_n).
	we also include len(t)/2+1 diffrent uniformly 
	spaced w's.
	"""
	p_min = np.max(t[3:]-t[:-3])
	w_max = np.pi*2/p_min
	p_max = np.max(t)-np.min(t)
	w_min = np.pi*2/p_max
	w_step = (w_max-w_min)*2/len(t)
	w = np.arange(w_min,w_max+w_step,w_step)
	return w	


def calcPeriodogram2(y,t,w,damp=0):
	"""Calculates the periodogram for the 
	time series data over a set of frequencies.
	To reduce problems with linear dependencies 
	(low rank) the frequencies are considered one
	at a time, and subsequent fits are performed on
	the residuals.  The periodogram is calculated 
	twice: first, the frequencies are considered in
	the initial order to estimate their weight; then,
	the second run is done in descending order.
	y	array - data vector
	t	array - time points corresponding to y
	w	array - frequencies to consider
	damp 	float - a piori damppining coef, defualt = 0, no dampping
	returns:
	B	periodogram weights corresponding to w 
	b	matrix of the actual fit coefficients, columns correspond to w
		coefficients in row -0 constant, -1 cos(), -2 sin().
	"""
	w = np.array(w)
	m = len(w)
	estB,estb = calcPeriodogram(y,t,w,False,damp)
	# sort results in decending order
	index = np.argsort(estB)[::-1]
	# sort are freqs
	tmpW = w[index]
	tmpB,tmpb = calcPeriodogram(y,t,tmpW,True,damp)
	b = np.zeros((m,3))
	B = np.zeros(m)
	# put back in right order
	B[index] = tmpB
	b[index,:] = tmpb
	return B,b

def calcPeriodogram(y,t,w,compensate=False,damp=0):
	"""Robustly Calculates the periodogram for the 
	time series data over a set of frequencies.
	To reduce problems with linear dependencies 
	(low rank) the frequencies are considered one
	at a time, and subsequent fits are performed on
	the residuals if requested.  The frequencies are considered in
	the same order as passed.  Robust linear regression
	is used to fit the model.  Currently the only option
	is iteratively reweighed least squares using the 
	Tukey's Biweight M-estimator.   
	y	array - data vector
	t	array - time points corresponding to y
	w	array - frequencies to consider
	compensate	bool - if ture subsequent frequancies
			are fitted against residuals
			Note: use if frequancies are
			passed in a "sensable" order
	damp	a priori dampining coef, defulat = 0, no damppining
	returns:
	B	periodogram weights corresponding to w 
	b	matrix of the actual fit coefficients, columns correspond to w
		coefficients in row -0 constant, -1 cos(), -2 sin().
	""" 
	n = len(w)
	b = np.zeros((n,3))
	B = np.zeros(n)
	for i in range(n):
		tmpb,res=cycleCoefRR(y,t,w[i],damp)
		if compensate: y = res
		b[i,:] = tmpb
		B[i] = float(n)/4.*(tmpb[1]**2+tmpb[2]**2)

	return B,b

def cycleCoefRR(y,t,w,damp=0):
	"""preforms a robust linear regression 
	on the time series data described. The linear model
	fit is to a sum of sin(w*t_i) and cos(w*t_i) and constant term 
	at the given frequency.
	Robust Regression is performed via iteratively reweighed
	least squares using the Tukey Biweight M-Estimator  
	y	array - vector of data
	t	array - time points corresponding to y
	w	scalar float - frequency to fit to
	returns:
	b	array - coefficients 0-const 1-cos 2-sin
	resid	fitted residual 
	"""
	n = len(t)
	X = np.ones((n,3))
	if np.abs(damp) < 1E-22:
		X[:,1] = np.cos(w*t)
		X[:,2] = np.sin(w*t)
	else:
		X[:,1] = np.exp(-damp*t)*np.cos(w*t)
		X[:,2] = np.exp(-damp*t)*np.sin(w*t)

	# Fit That Data!!!
	rlm = sm.RLM(y,X,M=norms.TukeyBiweight()).fit(conv='coefs',tol=1e-5)
	b = rlm.params
	resid = rlm.resid
	return b,resid

	
def estSignal(b,w,t,damp=0):
	"""Estimate a signal based 
	on the coefficents the frequancies 
	and the time points
	"""
	y = np.zeros(len(t))
	try:
		if np.abs(damp)<1E-22:
			for i in range(len(w)):
				y = y + b[i,0]+b[i,1]*np.cos(t*w[i])+b[i,2]*np.sin(t*w[i])
		else:
			for i in range(len(w)):
				y = y + b[i,0]+np.exp(-damp*t)*(b[i,1]*np.cos(t*w[i])+b[i,2]*np.sin(t*w[i]))

	except:
		if np.abs(damp)<1E-22:
			# maybe it was a scalar, ie single frequancy
			y = y + b[0]+np.exp(-damp*t)*(b[1]*np.cos(t*w)+b[2]*np.sin(t*w))
		else:
			y = y + b[0]+np.exp(-damp*t)*(b[1]*np.cos(t*w)+b[2]*np.sin(t*w))
	return y

def getFitTraj(y,t,w,tFit = [],damp=0):
	"""Estimate the fit based on allowed
	frequancies, w, and given data, y, 
	for time, t.  Then get the trajectory 
	for the fit at, tFit.
	damp is used for a priori damppining coef
	"""
	if len(tFit)==0:tFit = t
	
	B,b = calcPeriodogram2(y,t,w,damp)
	yFit = estSignal(b,w,tFit,damp)
	return(yFit)

def calcRCOD(yPrime,tPrime,w,damp=0,method='median'):
	"""Estimates the fit based on the 
	allowed freq in w, using robust regsiossion 
	and a priori dampening (if specified).  Then
	caclualtes and returns a simple robust
	coefficent of determination, which is 
	just the median squared error of the fit
	compared to a median squared error of the 
	null model (median of data)
	"""
	y,t = removenan(yPrime.copy(),tPrime.copy())
	B,b = calcPeriodogram2(y,t,w,damp)
	yFit = estSignal(b,w,t,damp)
	if method=='median':
		yNull = np.median(y)
		cod = (_medianAbsError(y,yNull)-_medianAbsError(y,yFit))/_medianAbsError(y,yNull)
	elif method=='mean':
		yVar = np.var(y)
		cod = (yVar-_meanSquaredError(y,yFit))/yVar
	else:
		raise ValueError('method '+method+' not avalible.')

	return(cod)


def compConstSlope(yPrime,tPrime,w,damp=0,method='mean'):
	"""Estimates the fit based on the 
	allowed freq in w, using robust regsiossion 
	and a priori dampening (if specified).  
	Calculates linear model wrt time as control model
	using the same robust regression 
	reports (mse(control) - mse(oscillator))/mse(control)
	"""
	y,t = removenan(yPrime.copy(),tPrime.copy())
	B,b = calcPeriodogram2(y,t,w,damp)
	yFit = estSignal(b,w,t,damp)

	X = np.ones((len(t),2))
	X[:,1] = t
	rlm = sm.RLM(y,X,M=norms.TukeyBiweight()).fit(conv='coefs',tol=1e-5)
	resid = rlm.resid
	if method == 'mean':
		eOsc = _meanSquaredError(y,yFit)	
		eCon = np.mean(resid**2)
	elif method=='median':
		eOsc = _medianAbsError(y,yFit)	
		eCon = np.median(np.abs(resid))
	else:
		raise ValueError('method '+method+' not avalible.')

	return((eCon-eOsc)/eCon)

def creatFreqVector(minP,maxP,n):
	"""given the mimimum and maximum period
	and the number of frequencies, create 
	a uniformly spaced frequency vector
	"""
	minW = np.pi*2/maxP
	maxW = np.pi*2/minP
	step = (maxW-minW)/(n-1.)
	w = np.arange(minW,maxW+step,step)
	return(w[:n])
	

	
def _medianAbsError(y,yHat):
	res = y - yHat
	return(np.median(np.abs(res)))
	
def _meanSquaredError(y,yHat):
	res = y - yHat
	return(np.mean(res**2))

def removenan(yPrime,tPrime):
	# this will remove nan values
	# since it has come up once and 
	# requiers relativly little work
	# I am includeing it here as a 
	# standard preprocessing step
	z = np.isnan(yPrime)
	y = np.array([])
	t = np.array([])
	for i in range(len(z)):
		if not z[i]:
			y = np.append(y,yPrime[i])
			t = np.append(t,tPrime[i])
		
	return y,t	
