import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np

def calcHPD(data, level=0.95):
	assert (0 < level < 1)	
	d = list(data)
	d.sort()	
	nData = len(data)
	nIn = int(round(level * nData))
	if nIn < 2 :
		raise RuntimeError("not enough data")	
	i = 0
	r = d[i+nIn-1] - d[i]
	for k in range(len(d) - (nIn - 1)) :
		rk = d[k+nIn-1] - d[k]
		if rk < r :
			r = rk
			i = k
	assert 0 <= i <= i+nIn-1 < len(d)	
	return [d[i], d[i+nIn-1]]
	

def hist(x,nbins = 50,xlab="x",ylab="Probability",fig_n= -1):
	# the histogram of the data
	if fig_n== -1: fig_n = int(np.random.uniform(0,10000))
	plt.figure(fig_n)
	n, bins, patches = plt.hist(x, nbins, normed=1, facecolor='green', alpha=0.5)
	# add a 'best fit' line
	# y = mlab.normpdf(bins, mu, sigma)
	# plt.plot(bins, y, 'r--')
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	#plt.title(r'histogram')
	# Tweak spacing to prevent clipping of ylabel
	plt.subplots_adjust(left=0.15)
	plt.show()
	

def plotXY(x,y,xlab="x",ylab="y",type="p",fig_n= -1):
	# the histogram of the data
	if fig_n== -1: fig_num = int(np.random.uniform(0,10000))
	else: fig_num = fig_n
	plt.figure(fig_num)
	if type == "p": p_type ='bo'
	if type == "l": p_type ='-'		
	p = plt.plot(x, y,p_type)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	#plt.title(r'histogram')
	# Tweak spacing to prevent clipping of ylabel
	plt.show()
	if fig_n == -1: plt.close()



def integrate_pdf(func,args,lower_lim=0,upper_lim=100,nbins=10000):
	# funtion, list of arguments
	v = np.linspace(lower_lim,upper_lim,nbins)
	d = v[1]-v[0]
	if upper_lim==0: return 0
	else: return sum(func(v,args))*d




















