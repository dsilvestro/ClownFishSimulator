import numpy as np
from numpy import *

def log_dgamma(L,a,b):  
	return ((a-1)*log(L)+(-b*L)-(log(b)*(-a)+ log(gamma(a))))
def log_dnorm(L,sd): 
	return (-(x**2/(2*sd**2)) - log(sd*sqrt(2*np.pi)))
def log_dcauchy(x,s):
	return (-log(np.pi*s * (1+ (x/s)**2)))

def dgamma(x,a,b):
	return( (1./b)**a * 1./gamma(a) * x**(a-1) *exp(-(1./b)*x))

def log_dunif(x,m,M): 
	if x>m and x<M: return( -log((M-m)) )
	else: return(-inf)


def update_parameter(i, m=-np.inf, M=np.inf, d=0.5): 
	ii = i+(np.random.random()-.5)*d
	if ii<m: ii=(ii-m)+m
	if ii>M: ii=(M-(ii-M))
	if ii<m: ii=i
	return ii

def update_parameter_normal(i, m=-np.inf, M=np.inf, d=0.5):
	ii = np.random.normal(i,d)
	if ii<m: ii=(ii-m)+m
	if ii>M: ii=(M-(ii-M))
	if ii<m: ii=i
	return ii

def update_multiplier_proposal(i,d):
	S=np.shape(i)
	u = np.random.uniform(0,1,S)
	l = 2*log(d)
	m = exp(l*(u-.5))
	ii = i * m
	return ii, sum(log(m))

