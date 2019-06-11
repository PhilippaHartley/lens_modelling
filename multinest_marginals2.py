#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
__doc__ = """
Script that does default visualizations (marginal plots, 1-d and 2-d).
Author: Johannes Buchner (C) 2013
"""
import numpy
from numpy import exp, log
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest
import matplotlib
matplotlib.rcParams.update({'font.size': 40})
matplotlib.rcParams['text.usetex'] 

if len(sys.argv) != 2:
	sys.stderr.write("""SYNOPSIS: %s <output-root> 
	output-root: 	Where the output of a MultiNest run has been written to. 
	            	Example: chains/1-
%s""" % (sys.argv[0], __doc__))
	sys.exit(1)

prefix = sys.argv[1]
print('model "%s"' % prefix)
if not os.path.exists(prefix + 'params.json'):
	sys.stderr.write("""Expected the file %sparams.json with the parameter names.
For example, for a three-dimensional problem:
["Redshift $z$", "my parameter 2", "A"]
%s""" % (sys.argv[0], __doc__))
	sys.exit(2)
parameters = json.load(open(prefix + 'params.json'))
n_params = len(parameters)

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s = a.get_stats()
m = numpy.array(s['marginals'][1]['5sigma'])
print( m)
s['marginals'][0]['5sigma'] = ((numpy.array(s['marginals'][0]['5sigma'])*2-964)*-1).tolist()
s['marginals'][1]['5sigma'] = (numpy.array(s['marginals'][1]['5sigma'])*2  -514).tolist()
s['marginals'][2]['5sigma'] = ((numpy.array(s['marginals'][2]['5sigma'])*2-964)*-1).tolist()
s['marginals'][3]['5sigma'] = (numpy.array(s['marginals'][3]['5sigma'])*2  -514).tolist()
s['marginals'][4]['5sigma'] = ((numpy.array(s['marginals'][4]['5sigma'])*2-964)*-1).tolist()
s['marginals'][5]['5sigma'] = (numpy.array(s['marginals'][5]['5sigma'])*2  -514).tolist()
s['marginals'][6]['5sigma'] = (numpy.array(s['marginals'][6]['5sigma'])*2).tolist()

s['marginals'][0]['3sigma'] = ((numpy.array(s['marginals'][0]['3sigma'])*2-964)*-1).tolist()
s['marginals'][1]['3sigma'] = (numpy.array(s['marginals'][1]['3sigma'])*2  -514).tolist()
s['marginals'][2]['3sigma'] = ((numpy.array(s['marginals'][2]['3sigma'])*2-964)*-1).tolist()
s['marginals'][3]['3sigma'] = (numpy.array(s['marginals'][3]['3sigma'])*2  -514).tolist()
s['marginals'][4]['3sigma'] = ((numpy.array(s['marginals'][4]['3sigma'])*2-964)*-1).tolist()
s['marginals'][5]['3sigma'] = (numpy.array(s['marginals'][5]['3sigma'])*2  -514).tolist()
s['marginals'][6]['3sigma'] = (numpy.array(s['marginals'][6]['3sigma'])*2).tolist()

s['marginals'][0]['1sigma'] = ((numpy.array(s['marginals'][0]['1sigma'])*2-964)*-1).tolist()
s['marginals'][1]['1sigma'] = (numpy.array(s['marginals'][1]['1sigma'])*2  -514).tolist()
s['marginals'][2]['1sigma'] = ((numpy.array(s['marginals'][2]['1sigma'])*2-964)*-1).tolist()
s['marginals'][3]['1sigma'] = (numpy.array(s['marginals'][3]['1sigma'])*2  -514).tolist()
s['marginals'][4]['1sigma'] = ((numpy.array(s['marginals'][4]['1sigma'])*2-964)*-1).tolist()
s['marginals'][5]['1sigma'] = (numpy.array(s['marginals'][5]['1sigma'])*2  -514).tolist()
s['marginals'][6]['1sigma'] = (numpy.array(s['marginals'][6]['1sigma'])*2).tolist()

s['marginals'][0]['median'] = ((numpy.array(s['marginals'][0]['median'])*2-964)*-1).tolist()
s['marginals'][1]['median'] = (numpy.array(s['marginals'][1]['median'])*2  -514).tolist()
s['marginals'][2]['median'] = ((numpy.array(s['marginals'][2]['median'])*2-964)*-1).tolist()
s['marginals'][3]['median'] = (numpy.array(s['marginals'][3]['median'])*2  -514).tolist()
s['marginals'][4]['median'] = ((numpy.array(s['marginals'][4]['median'])*2-964)*-1).tolist()
s['marginals'][5]['median'] = (numpy.array(s['marginals'][5]['median'])*2  -514).tolist()
s['marginals'][6]['median'] = (numpy.array(s['marginals'][6]['median'])*2).tolist()



#= (numpy.array(s['modes']['sigma'][0])*2-964).tolist()
s['modes'][0]['sigma'][0] = s['modes'][0]['sigma'][0]*2
s['modes'][0]['sigma'][1] = s['modes'][0]['sigma'][1]*2
s['modes'][0]['sigma'][2] = s['modes'][0]['sigma'][2]*2
s['modes'][0]['sigma'][3] = s['modes'][0]['sigma'][3]*2
s['modes'][0]['sigma'][4] = s['modes'][0]['sigma'][4]*2
s['modes'][0]['sigma'][5] = s['modes'][0]['sigma'][5]*2
s['modes'][0]['sigma'][6] = s['modes'][0]['sigma'][6]*2




s['modes'][0]['mean'][0] = (s['modes'][0]['mean'][0]*2-964)*-1
s['modes'][0]['mean'][1] = s['modes'][0]['mean'][1]*2-514
s['modes'][0]['mean'][2] = (s['modes'][0]['mean'][2]*2-964)*-1
s['modes'][0]['mean'][3] = s['modes'][0]['mean'][3]*2-514
s['modes'][0]['mean'][4] = (s['modes'][0]['mean'][4]*2-964)*-1
s['modes'][0]['mean'][5] = s['modes'][0]['mean'][5]*2-514
s['modes'][0]['mean'][6] = s['modes'][0]['mean'][6]*2


print('modes ', s['modes'][0]['sigma'])








#s['modes'][1]['sigma'] = (numpy.array(s['modes'][1]['sigma'])*2  -514).tolist()

#s['modes'][0]['mean'] = (numpy.array(s['modes'][0]['mean'])*2 -964).tolist()
#s['modes'][1]['mean'] = (numpy.array(s['modes'][1]['mean'])*2  -514).tolist()



json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)


print('  marginal likelihood:')
print('    ln Z = %.1f +- %.1f' % (s['global evidence'], s['global evidence error']))
print('  parameters:')
for p, m in zip(parameters, s['marginals']):
	lo, hi = m['1sigma']
	med = m['median']
	print (lo, hi, med)
	sigma = numpy.abs(hi - lo) / 2
	print (sigma)
	i = max(0, int(-numpy.floor(numpy.log10(sigma))) + 1)
	fmt = '%%.%df' % i
	fmts = '\t'.join(['    %-15s' + fmt + " +- " + fmt])
	print(fmts % (p, med, sigma))

print('creating marginal plot ...')
p = pymultinest.PlotMarginal(a)

values = a.get_equal_weighted_posterior()
print( numpy.std(values[:,4]))
values[:,0]=values[:,0]*2-964
values[:,1]=values[:,1]*2-514
values[:,2]=values[:,2]*2-964
values[:,3]=values[:,3]*2-514
values[:,4]=values[:,4]*2-964
values[:,5]=values[:,5]*2-514
values[:,6]=values[:,6]*2

values[:,0]*=-1
values[:,2]*=-1
values[:,4]*=-1


assert n_params == len(s['marginals'])
modes = s['modes']

dim2 = os.environ.get('D', '1' if n_params > 20 else '2') == '2'
nbins = 100 if n_params < 3 else 20
if dim2:
	plt.figure(figsize=(7*n_params, 7*n_params))
	plt.rc('text', usetex=True)
	for i in range(n_params):
		cch = 'crimson'
		ax = plt.subplot(n_params, n_params,(n_params*i)+1+i)
		ax.set_xlabel(parameters[i])
		lower = modes[0]['mean'][i] - 3*modes[0]['sigma'][i] 
		upper = modes[0]['mean'][i] + 3*modes[0]['sigma'][i]	
		m = s['marginals'][i]
		print ('m: ')
		print (m)
		ax.set_xlim(lower, upper)
		
		print ('u, l ', upper , lower)
		for mo in modes:
			print (mo['mean'][i] - 3*mo['sigma'][i] )
			print( mo['mean'][i] + 3*mo['sigma'][i] )
		print( 'xlim: ')
		print(m['5sigma'])
		oldax= plt.gca()
		print( values[:,i])
		x,w,patches = oldax.hist(values[:,i], bins=nbins, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
		oldax.set_ylim(0, x.max())
	
		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
	#	p.plot_marginal(i, ls='-', color=cch, linewidth=4)
		newax.set_ylim(0, 1)
	
		ylim = newax.get_ylim()
		print( 'ylim: ')
		print(ylim)
		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
	#	center = m['mean']
		center = s['modes'][0]['mean'][i] 
		low1, high1 = m['1sigma']
		print('center ', center, low1, high1)
		newax.errorbar(x=center, y=y,
			xerr=numpy.transpose([[center - low1, high1 - center]]), 
			color=cch, linewidth=3, marker='s')
		oldax.set_yticks([])
        

		#oldax.set_xticklabels([])
		#newax.set_yticks([])
		plt.rc('text', usetex=True)
		if i == 0:
		
		    plt.ylabel("SX1")
		else:
			#pass
		    newax.set_yticklabels([])      
		if not i == n_params-1:
			newax.set_xlabel("GPA")
			newax.set_xticklabels([])    
		lab = newax.get_xticklabels()
		rot = 60
		for k in lab:
			k.set_rotation(rot) 
		lab2 = oldax.get_xticklabels()
		for k in lab2:
			k.set_rotation(rot) 	
		ylim = oldax.get_ylim()
		lower = modes[0]['mean'][i] - 3*modes[0]['sigma'][i] 
		upper = modes[0]['mean'][i] + 3*modes[0]['sigma'][i] 
		newax.set_xlim(lower, upper)
		print( 'xlim: ')
		print(m['3sigma'])
		si = m['5sigma']
		oldax.set_xlim(lower, upper)
		#plt.close()
		plt.xticks(rotation=rot)	
		for j in range(i):
			sp = n_params*(i)+1+j
			ax = plt.subplot(n_params, n_params, (n_params*i)+j+1)
			#ax.xaxis.set_tick_params(rotation=45)
	#		if not isinstance((n_params * (j + 1) + i + 1)/n_params-j, (int, numpy.integer)):
#				ax.set_xticklabels([])#
		#		ax.set_yticklabels([])
		#	if i>j+1:
	#			ax.set_xticklabels([])
	#			ax.set_yticklabels([])
				
			#ax.set_aspect('equal')
			#ax.set_xlim(si)
			lower = modes[0]['mean'][j] - 3*modes[0]['sigma'][j] 
			upper = modes[0]['mean'][j] + 3*modes[0]['sigma'][j]	
			print('u, l:, ', lower, upper)
			if i == n_params-1:
			
				p.plot_conditional(j, i, bins=20, cmap = plt.cm.gray_r)
			else:
				p.plot_conditional(j, i, bins=20, cmap = plt.cm.gray_r)
			for mo in modes:
				print ('from modes ' , mo['mean'][j],mo['mean'][i], mo['sigma'][j])
				plt.errorbar(x=mo['mean'][j], y=mo['mean'][i], xerr=mo['sigma'][j], yerr=mo['sigma'][i], color = cch, linewidth = 4)
				print ('u, l')
				print (mo['mean'][j] - 3*mo['sigma'][j] )
				print( mo['mean'][j] + 3*mo['sigma'][j] )
			plt.xlim(lower, upper)	
				
			xl = plt.xlim()
			print ('xl: ', xl)
			#sys.exit()	
			#if i==j+1:
			if j==0:
				plt.ylabel(parameters[i])
			else: 
			   # pass
			    ax.set_yticklabels([])
			plt.xticks(rotation=rot)	
			plt.xlabel(parameters[j])	
			print (parameters[j])
			#plt.savefig('cond_%s_%s.pdf' % (parameters[i], parameters[j]), bbox_tight=True)
		#	plt.close()
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig(prefix + 'marg.pdf')
	plt.savefig(prefix + 'marg.png')
	plt.close()
else:
	from matplotlib.backends.backend_pdf import PdfPages
	sys.stderr.write('1dimensional only. Set the D environment variable \n')
	sys.stderr.write('to D=2 to force 2d marginal plots.\n')
	pp = PdfPages(prefix + 'marg1d.pdf')
	
	for i in range(n_params):
		plt.figure(figsize=(3, 3))
		plt.xlabel(parameters[i])
		plt.locator_params(nbins=5)
		
		m = s['marginals'][i]
		iqr = m['q99%'] - m['q01%']
		xlim = m['q01%'] - 0.3 * iqr, m['q99%'] + 0.3 * iqr
		#xlim = m['5sigma']
		plt.xlim(xlim)
	
		oldax = plt.gca()
		x,w,patches = oldax.hist(values[:,i], bins=numpy.linspace(xlim[0], xlim[1], 20), edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
		oldax.set_ylim(0, x.max())
	
		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
		newax.set_ylim(0, 1)
	
		ylim = newax.get_ylim()
		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
		center = m['median']
		low1, high1 = m['1sigma']
		#print center, low1, high1
		newax.errorbar(x=center, y=y,
			xerr=numpy.transpose([[center - low1, high1 - center]]), 
			color='blue', linewidth=2, marker='s')
		oldax.set_yticks([])
		newax.set_ylabel("Probability")
		ylim = oldax.get_ylim()
		newax.set_xlim(xlim)
		oldax.set_xlim(xlim)
		plt.savefig(pp, format='pdf', bbox_inches='tight')
		plt.close()
	pp.close()


