'''---------------------------------------------------------------------------------------
HiCPlotter: plotting Hi-C data with additional datasets
------------------------------------------------------------------------------------------'''

import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " HiCPlotter needs python2.7.*!\n"
	sys.exit()

import platform
if platform.platform().split('-')[0]=='Linux' or platform.platform().split('-')[0]=='Windows':
	import matplotlib
	matplotlib.use('Agg')

from math  import sqrt, isnan, floor, ceil, pi
from numpy import log2, array, max
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker  import MultipleLocator
from matplotlib.patches import Polygon, Rectangle, Circle
from scipy.signal import argrelextrema
from scipy import ndimage
import scipy.sparse as sps
import matplotlib.pyplot as plt
import numpy as np
import argparse
import bisect
import warnings
import logging

version = "0.7.0"

def read_HiCdata(filename,header=0,footer=0,clean_nans=True,smooth_noise=0.5,ins_window=5,rel_window=8,plotInsulation=True,plotTadDomains=False,randomBins=False):
	
	'''
    load Hi-C interaction matrix from text file
    
    parameters:

    filename: file name. format "chr\tstart\tend\tdata1\tdata2\t..."
    clean_nans: replace nan with 0's in rows/columns from all matrix
    smooth_noise: variable values under will be replace with 0's to clean noise in the matrix
    ins_window: window size for scanning diagonal of the matrix for insulation scores
    rel_window: relative extrama extension size - will be extend to both directions
    
    returns:
    
    matrix: data matrix over the selected set of chromosome.
    nums: insulation scores array.
    tricks: putative insulator sites

    '''
	
	try:
		matrix = np.genfromtxt(filename,skip_header=header,filling_values='0',skip_footer=footer) 
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	if not randomBins and len(matrix[:,1]) != len(matrix[1,:]):
		print len(matrix[:,1]),len(matrix[1,:])
		print >>sys.stderr, 'unbalanced matrix('+filename+')! input should be a square matrix'
		raise SystemExit
	
	if plotInsulation or plotTadDomains and not randomBins: nums,tricks=insulation(matrix,ins_window,rel_window)
	else: nums=[];tricks=[];
	
	if clean_nans: matrix[np.isnan(matrix)]=0
	matrix[matrix<smooth_noise]=0
	return matrix,nums,tricks


def read_sparseHiCdata(filename,chromosome,bedFile,startBin,endBin,wholeGenome=False,smooth_noise=0.5,ins_window=5,rel_window=8,plotInsulation=True,plotTadDomains=False,randomBins=False):
						
	'''
    load Hi-C interaction matrix from triple-column sparse file
    
    parameters:

    filename: file name. format "bin1\tbin2\tdata\n..."
    chromosome: plotting which chromosome
    bedFile: a bed file for locations of bins
    startBin: starting bin - 0 zero-based
    endBin: end point for the plot
    wholeGenome: for plotting more than one chromosome interactions. chromosome parameter will be used for until which chromosome interactions will be plotted.
    smooth_noise: variable values under will be replace with 0's to clean noise in the matrix
    ins_window: window size for scanning diagonal of the matrix for insulation scores
    rel_window: relative extrama extension size - will be extend to both directions
    
    returns:
    
    matrix: data matrix over the selected set of chromosome.
    nums: insulation scores array.
    tricks: putative insulator sites

    '''
	
	
	chromosomes = {}
	
	try:
		bed = open(bedFile,'r') 
	except IOError:
		print >>sys.stderr, 'cannot open', bedFile
		raise SystemExit
	
	for line in bed.readlines():
		tags = line.strip().split("\t")
		if tags[0]=='chrM':continue
		if tags[0] not in chromosomes.keys():
			chromosomes[tags[0]]=[]
			chromosomes[tags[0]].append(int(tags[3]))
		else: chromosomes[tags[0]].append(int(tags[3]))
	
	if not wholeGenome:
		clast = chromosomes[chromosome][-1]-1
		start = chromosomes[chromosome][0]+startBin
		end = chromosomes[chromosome][0]+endBin
	
		end=clast if end == chromosomes[chromosome][0] else end
		if end > clast: end=clast
		if start > clast: start=chromosomes[chromosome][0]
	else:
		start = 1
		end = chromosomes[chromosome][-1]
		clast = end
	
	length = end-start+1
	
	mtx = sps.dok_matrix((length, length), dtype=np.int)
	
	try:
		matrixFile = open(filename,'r') 
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	for line in matrixFile.xreadlines():
		tags = line.strip().split("\t")
		if int(tags[0]) <= end and int(tags[0])>=start :
			if int(tags[1]) <= end and int(tags[1])>=start :
				mtx[int(tags[0])-start, int(tags[1])-start] = int(round(float(tags[2])))
				mtx[int(tags[1])-start, int(tags[0])-start] = int(round(float(tags[2])))
		if int(tags[0]) > end: break
	
	matrix = mtx.todense()
	
	if plotInsulation or plotTadDomains and not wholeGenome: nums,tricks=insulation(matrix,ins_window,rel_window,True,startBin)
	else: nums=[];tricks=[];
	
	#matrix[matrix<smooth_noise]=0
	return matrix,nums,tricks,clast-chromosomes[chromosome][0]+1

def read_bedGraph(filename,resolution,chromosome): # add stopping after certain chromosome passed
	
	'''
    reads bedGraph files for various file type plottings

    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	
	returns:
	x_scores = location along the given chromosome - start sites
	x_scores2 = location along the given chromosome - end sites
	y_scores = signal scores for the assay
	colors = allow for colors option
    '''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	x_scores=[]
	x_scores2=[]
	y_scores=[] 
	colors=[]
	texts=[]
	
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome:
			x_scores.append(float(tags[1])/resolution)
			x_scores2.append(float(tags[2])/resolution)
			if len(tags) > 3:
				y_scores.append(float(tags[3]))
				if len(tags) > 4:
					hex = '#%02x%02x%02x' % (int(tags[4].split(',')[0]), int(tags[4].split(',')[1]), int(tags[4].split(',')[2]))
					colors.append(hex)
					if len(tags) > 5:
						texts.append(tags[5])
				
	if len(y_scores) !=0 and len(y_scores)!=len(x_scores):
		print >>sys.stderr, 'BedGraph('+filename+') has some missing values'
		raise SystemExit
	if len(x_scores)==0 or len(x_scores2)==0:
		print >>sys.stderr, 'BedGraph('+filename+') has some missing values'
		raise SystemExit
	# color and text controls
	return x_scores,x_scores2,y_scores,colors,texts
	
def read_peakFile(filename,resolution,chromosome): # add stopping after certain chromosome passed
	
	'''
    reads peak files for annotating the matrix

    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	
	returns:
	origin_x = location along x axis on the given chromosome
	origin_y = location along y axis on the given chromosome
	radius = radius of circle
	colors = allow for colors option
    '''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	origin_x=[]
	origin_y=[]
	radius=[] 
	colors=[]
	
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome and tags[3]==chromosome:
			x1 = float(tags[1])/resolution
			x2 = float(tags[2])/resolution
			origin_x.append(x1+(x2-x1)/2)
			radius.append((x2-x1)/2)
			y1 = float(tags[4])/resolution
			y2 = float(tags[5])/resolution
			origin_y.append(y1+(y2-y1)/2)
			
			if len(tags) > 5:
				hex = '#%02x%02x%02x' % (int(tags[6].split(',')[0]), int(tags[6].split(',')[1]), int(tags[6].split(',')[2]))
				colors.append(hex)
			
	if len(origin_y) !=0 and len(origin_x)!=len(origin_y):
		print >>sys.stderr, 'Peak file ('+filename+') has some missing values'
		raise SystemExit
	if len(origin_x)==0 or len(origin_y)==0:
		print >>sys.stderr, 'Peak file ('+filename+') has some missing values'
		raise SystemExit
	# color control
	return origin_x,origin_y,radius,colors

def read_epilogos(filename,resolution,chromosome,start,end): # add stopping after certain chromosome passed
	
	'''
    reads epilogos file format (http://compbio.mit.edu/epilogos/)
	# you can download epilogos files for human genome from http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/epilogos/
    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	start: starting bin - 0 zero-based
    end: end point for the plot
	
	returns:
	x_scores = location along the given chromosome - start sites
	y_dict = a dictionary containing enrichments of each states for a given location
    '''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	x_scores=[]
	y_dict={} 
	start = resolution * start
	end = resolution * end
	
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome and int(tags[1])>start and int(tags[2]) < end:
			x_scores.append(float(tags[1])/resolution)
			cols = np.array(tags[3].split(':')[2].split(','))
			for x in range(1,len(cols)):
				if x % 2 == 1 and cols[x].replace('[','').replace(']','').replace(' ','') not in y_dict.keys():
					y_dict[cols[x].replace('[','').replace(']','').replace(' ','')]=[]
					y_dict[cols[x].replace('[','').replace(']','').replace(' ','')].append(float(cols[x-1].replace('[','').replace(']','').replace(' ','')))
				elif x % 2 == 1 :
					y_dict[cols[x].replace('[','').replace(']','').replace(' ','')].append(float(cols[x-1].replace('[','').replace(']','').replace(' ','')))
					
	if len(y_dict.keys()) ==0:
		print >>sys.stderr, 'Epilogos File('+filename+') has some missing values'
		raise SystemExit
	return x_scores,y_dict

def read_genes(filename,resolution,chromosome,start,end):
	
	'''
    reads a sorted gene location file

    parameters:

    filename: file name. format could be either "chr\tstart\tend" or "chr\tstart\tend\tvalue..."
    resolution: bin size for the matrix
	start: starting bin - 0 zero-based
    end: end point for the plot
    
	returns:
	genes : a dictionary for each gene and respective plot locations
	row_count: a number to indicate how many rows will be used
	row_genes: a dictionary for end locations of the genes on each row
	'''
	
	try:
		fone=open(filename,'r')
	except IOError:
		print >>sys.stderr, 'cannot open', filename
		raise SystemExit
	
	start = resolution * start
	end = resolution * end
	
	minDist = 8000
	
	genes = {}
	row_list = []
	row_genes = {}
	current_start=0;current_end=0;prev_end=0;
	for line in fone.xreadlines():
		tags = line.strip().split("\t")
		if tags[0]==chromosome:
			if int(tags[1]) >= start and int(tags[2]) <= end and tags[1]+'-'+tags[2] not in genes.keys():
				
				if len(row_list)==0:
					current_start = int(tags[1])
					current_end = int(tags[2])
					prev_start = int(tags[1])
					genes[tags[1]+'-'+tags[2]]=[]
					genes[tags[1]+'-'+tags[2]].append(1)
					genes[tags[1]+'-'+tags[2]].append(tags[3])
					row_list.append(current_end)
					row_genes[1]=[]
					row_genes[1].append(current_end)
				else:
					if prev_start > int(tags[1]):
						print prev_end, int(tags[1])
						print >>sys.stderr, 'Gene File ('+filename+') is not sorted.'
						raise SystemExit
					else:
						current_end = int(tags[2])
						current_start = int(tags[1])
						execute=0
						genes[tags[1]+'-'+tags[2]]=[]
						for item in range(0,len(row_list)):
							if current_start > row_list[item]+minDist: 
								row_list[item]=current_end
								execute=1
								genes[tags[1]+'-'+tags[2]].append(item+1)
								if item+1 not in row_genes.keys(): row_genes[item+1]=[]
								row_genes[item+1].append(current_end)
								break
						if execute == 0:
							genes[tags[1]+'-'+tags[2]].append(len(row_list)+1)
							row_list.append(current_end)
							if len(row_list) not in row_genes.keys(): row_genes[len(row_list)]=[]
							row_genes[len(row_list)].append(current_end)
							
						genes[tags[1]+'-'+tags[2]].append(tags[3])
						if len(tags)>5:
							genes[tags[1]+'-'+tags[2]].append(tags[4])
							genes[tags[1]+'-'+tags[2]].append(tags[5])
							genes[tags[1]+'-'+tags[2]].append(tags[6])
							if len(tags)>7:
								hex = '#%02x%02x%02x' % (int(tags[7].split(',')[0]), int(tags[7].split(',')[1]), int(tags[7].split(',')[2]))
								genes[tags[1]+'-'+tags[2]].append(hex)
								
					prev_start = current_start
							
	if len(genes.keys()) ==0:
		print >>sys.stderr, 'Gene File ('+filename+') has some missing values'
		raise SystemExit
	
	return genes,len(row_list)+1,row_genes

def where(start,end,arr):
    """Find where the start location and end location indexes in an array"""
    
    astart = bisect.bisect_left(arr, start)
    aend = bisect.bisect_right(arr[start:], end) + start
        
    return astart, aend

def get_ellipse_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2):
    """ Draws an ellipse using (360*k + 1) discrete points
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees clockwise from positive x-axis
        
    this function is obtained from : http://scipy-central.org/item/23/2/plot-an-ellipse
    """
    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])
 
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    
    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts


def insulation(matrix,w=5,tadRange=10,triple=False,mstart=0):
	
	'''
    calculate relative minima in a given matrix

    parameters:

    matrix: data matrix over the selected set of chromosomes.
    start: retain after x-th bin.
    end: continues until x-th bin.
	w: window size for scanning diagonal of the matrix for insulation scores
    tadRange: relative extrama extension size - will be extend to both directions
    
	returns:
	nums: insulation scores array.
    tricks: putative insulator sites
	
    '''
    
	start=0;end=len(matrix)-1
	scores = []
	indexes = []
	pBorders=[]
	
	
	for i in xrange(start,end,1):
		diag=0;counter=0
		for j in xrange(i,i+w):
			if j == end: break
			else:
				if triple:
					if isnan(matrix[j,j-2*counter]): diag +=0 # pad with zeros for nan
					else: diag += matrix[j,j-2*counter]
				else:
					if isnan(matrix[j][j-2*counter]): diag +=0 # pad with zeros for nan
					else: diag += matrix[j][j-2*counter]
				counter+=1
		scores.append(diag)
		indexes.append(i)

	arr= np.array(scores)

	arr[arr == 0] = 10000
	borders = argrelextrema(arr, np.less,order=tadRange) #also try this function scipy.signal.find_peaks_cwt
	regions = borders[0]
	
	
	for item in range(0,len(regions)-1):	
		if item == 0: #check the first boundary
			if len(np.nonzero(scores[slice(regions[item],regions[item+1])])[0]) < regions[item+1]-regions[item]-w/2-1:
				current = np.nonzero(scores[slice(0,regions[item])])[0]
				if current[0] not in pBorders and current[0]+1 not in regions: pBorders.append(current[0])
				if regions[item] not in pBorders and regions[item]+1 not in regions: pBorders.append(regions[item])
				current = np.nonzero(scores[slice(regions[item],regions[item+1])])[0]
				if regions[item]+current[-1] not in pBorders and regions[item]+current[-1]+1 not in regions: pBorders.append(regions[item]+current[-1])
				elif regions[item]+current[-1] not in pBorders and regions[item]+current[-1]+1 in regions: pBorders.append(regions[item]+current[-2]) # to get closer bin but needs to be improved
			else:
				if regions[item] not in pBorders and regions[item]+1 not in regions: pBorders.append(regions[item])
		else:
			current = np.nonzero(scores[slice(regions[item],regions[item+1])])[0]
			if len(current) >= regions[item+1]-regions[item]-w/2-1:
				if regions[item] not in pBorders and regions[item]+1 not in regions: pBorders.append(regions[item])
			elif len(current) < regions[item+1]-regions[item]-w/2-1:
				current = np.nonzero(scores[slice(regions[item],regions[item]+tadRange)])[0]
				if regions[item] not in pBorders and regions[item]+1 not in regions: pBorders.append(regions[item])
				if regions[item]+current[0] not in pBorders and regions[item]+current[0]+1 not in regions: pBorders.append(regions[item]+current[0])
				current = np.nonzero(scores[slice((regions[item]+tadRange),regions[item+1])])[0]
				if len(current)==0:
					if regions[item+1] not in pBorders and regions[item]+1 not in regions: pBorders.append(regions[item+1])
				else:
					if regions[item]+tadRange+current[0] not in pBorders and regions[item]+tadRange+current[0]+1 not in regions: pBorders.append(regions[item]+tadRange+current[0])
					if regions[item+1] not in pBorders and regions[item+1]+1 not in regions: pBorders.append(regions[item+1])
	if regions[-1] not in pBorders: pBorders.append(regions[-1])
	if len(matrix) not in pBorders: pBorders.append(len(matrix))
	
	if triple: pBorders=map(lambda x:x+mstart, pBorders)
	
	return scores, pBorders

def HiCplotter(files=[],names=[],resolution=100000,chromosome='',output='',histograms=[],histLabels=[],fillHist=[],histMax=[],verbose=False,fileHeader=0,fileFooter=0,matrixMax=0,histColors=[],barPlots=[],barLabels=[],plotGenes='',superImpose=False,\
			start=0,end=0,tileLabels=[],tilePlots=[],tileColors=[],tileText=False,arcLabels=[],arcPlots=[],arcColors=[],peakFiles=[],epiLogos='',window=5,tadRange=8,tripleColumn=False,bedFile='',barColors=[],dPixels=200,compareEx='',compareSm='',\
			smoothNoise=0.5,cleanNANs=True,plotTriangular=True,plotTadDomains=False,randomBins=False,wholeGenome=False,plotPublishedTadDomains=False,plotDomainsAsBars=False,imputed=False,barMax=[],spine=False,plotDomainTicks=True,triangularHeight=False,\
			highlights=0,highFile='',heatmapColor=3,highResolution=True,plotInsulation=True,plotCustomDomains=False,publishedTadDomainOrganism=True,customDomainsFile=[],compare=False,pair=False,domColors=[],oExtension='',geneLabels=True,dark=False):
	
	'''
    plot the interaction matrix with additional datasets

 	Required parameters:

    files 			(-f)		: a list of filenames to be plotted.
    name 			(-n) 		: a list of labels for the experiment.
    chr				(-chr)		: chromosome to be plotted.
    output			(-o)		: prefix for the output file.
    
    Optional parameters:
    
    verbose			(-v)		: print version and arguments into a file.
    tripleColumn	(-tri)		: a boolean if input file is from HiC-Pro pipeline.
    dark			(-d)		: a boolean to use black background for the output.
    bedFile			(-bed)		: a file name for bin annotations, if -tri parameter is set.
    plotGenes		(-g)		: a sorted bed file for plotting the locations of the genes.
    geneLabels		(-gl)		: a boolean for plotting gene labels (1:default) or not (0).
    histograms		(-hist)		: a list of filenames to be plotted as histogram.
    histLabels		(-h)		: a list of labels for the histograms.
    fillHist		(-fhist)	: a list whether each histogram will be filled (1) or not (0:default).
    histColors		(-hc)		: a list of hexadecimal number for histogram filling colors.
    histMax 		(-hm)		: a list of integer for maximum values of histograms.
    superImpose		(-si)		: a boolean to overlap two histogram files inside the same track (default:0) enable(1).
    start			(-s)		: retain after x-th bin (0:default).
    end				(-e)		: continues until x-th bin (default: length of the matrix).
    resolution		(-r)		: resolution of the bins (default: 100000).
    matrixMax		(-mm)		: an integer value for the interaction matrix heatmap scale upper-limit.
    barPlots		(-b)		: a list of filenames to be plotted as bar plots.
    barLabels		(-bl)		: a list of labels for the bar plots.
    barColors		(-bc)		: a list of hexadecimal numbers for coloring the bar plots.
    barMax	 		(-bm)		: a list of integer for maximum values of bar plots.
    tilePlots		(-t)		: a list of filenames to be plotted as tile plots.
    tileLabels		(-tl)		: a list of labels for the tile plots.
    tileColors		(-tc)		: a list of hexadecimal numbers for coloring the tile plots.
    tileText		(-tt)		: a boolean whether text will be displayed above tiles (0:default) or not (1).
    arcPlots		(-a)		: a list of filenames to be plotted as arc plots.
    arcLabels		(-al)		: a list of labels for the arc plots.
    arcColors		(-ac)		: a list of hexadecimal numbers for coloring the arc plots.
    highlights		(-high)		: a boolean for enabling highlights on the plot (0:default), enable(1). 
    highFile		(-hf)		: a file name for a bed file to highlight selected intervals.
    peakFiles 		(-peak)		: a list of filenames to be plotted on the matrix.
    epiLogos 		(-ep)		: a filename to be plotted as Epilogos format.
    oExtension 		(-ext)		: an extension name for the output file format - default jpeg.
    imputed 		(-im)		: a boolean if imputed epilogos will be plotted. (default:0 for observed)
    spine			(-spi)		: a boolean to remove top and left borders for each tracks (default:0) enable(1).
    compare			(-c)		: a boolean to plot log2 compare first two matrices (default:0) enable(1).
    compareExt		(-ce)		: comma separated two integers for log2 comparison matrix color spectrum (e.g. 2,4 for -2 to 4). 
    compareSm		(-cs)		: comma separated two integers for log2 comparison matrix smoothing (e.g. for 0,2 values between 0-2 will be white in image).
    pair			(-p)		: a boolean to plot log2 pair-wise matrix comparisons (default:0) enable(1).
    window			(-w)		: an integer of distance to calculate insulation score.
    tadRange		(-tr)		: an integer of window to calculate local minima for TAD calls.
    fileHeader		(-fh)		: an integer for how many lines should be ignored in the matrix file (0:default).
    fileFooter		(-ff)		: an integer for how many lines should be skipped at the end of the matrix file (0:default).
    smoothNoise		(-sn)		: a floating-point number to clean noise in the data.
    heatmapColor	(-hmc)		: an integer for choosing heatmap color codes: Greys(0), Reds(1), YellowToBlue(2), YellowToRed(3-default), Hot(4), BlueToRed(5).
    cleanNANs		(-cn)		: a boolean for replacing NaNs in the matrix with zeros (1:default) or not (0).
    plotTriangular	(-ptr)		: a boolean for plotting rotated half matrix (1:default) or not (0).
    domColors		(-dc)		: a list of hexadecimal numbers for coloring the domain plots.
    plotTadDomains	(-ptd)		: a boolean for plotting TADs identified by HiCPlotter (1) or not (0:default).
    plotPublishedTadDomins	(-pptd)	: a boolean for plotting TADs from Dixon et, al. 2012 (1:default) or not (0).
    plotDomainsAsBars		(-ptdb)	: a boolean for plotting TADs as bars (1) instead of triangles (0:default)
    highResolution	(-hR)		: a boolean whether plotting high resolution (1:default) or not (0).
    dPixels			(-dpi)		: an integer to determine dots per inch in matrix, higher values for higher resolution (default:200).
    plotInsulation	(-pi)		: a boolean for plotting insulation scores (0:default) or plot (1).
    randomBins		(-rb)		: a boolean for plotting random resolution data (1:default) or not (0).
    wholeGenome		(-wg)		: a boolean for plotting whole genome interactions (1:default) or not (0).
    plotCustomDomains		(-pcd)	: a list of file names to be plotted beneath the matrix.
    publishedTadDomainOrganism 	(-ptdo)	: a boolean for plotting human (1:default) or mouse (0) TADs from Dixon et, al. 2012.
    customDomainsFile			(-pcdf)	: a list of filenames to be plotted as TADs for each experiments.
	
	'''
	
	
	numOfcols = len(files)
	numOfrows = 4
	if plotTriangular: numOfrows+=2
	if len(plotGenes)>0: numOfrows+=2
	if plotTadDomains: numOfrows+=1
	if plotInsulation: numOfrows+=1
	if epiLogos: numOfrows+=1
	if len(histograms)>0 and not superImpose: numOfrows+=len(histograms[0].split(','))
	if superImpose : numOfrows+=1
	if len(barPlots)>0: numOfrows+=len(barPlots[0].split(','))
	if len(tilePlots)>0: numOfrows+=len(tilePlots[0].split(','))
	if len(arcPlots)>0: numOfrows+=len(arcPlots[0].split(','))
	if plotCustomDomains or plotPublishedTadDomains and not plotTadDomains: numOfrows+=1
	
	if compare and not pair: numOfcols+=1; files.append('pseudo'); matrix1=[];matrix2=[]
	if compare and pair: numOfrows+=(len(files)-1)*4
	if pair: marray=[]
	
	fig=plt.figure(figsize=(numOfcols*5+2.5, numOfrows+numOfrows/2+0.5), facecolor='w', edgecolor='w')
	fig.set_size_inches(numOfcols*5+2.5, numOfrows+numOfrows/2+0.5)
	if superImpose: fig.subplots_adjust(hspace=0.48,wspace=1.25)
	else: fig.subplots_adjust(hspace=0.48,wspace=1.0)
	if dark : plt.style.use('dark_background')
	
	ymaxlims = []
	yminlims = []
	cmatrix = 0
	ins_score = 0
	mlength = 0
	
	cmaps = ['Greys','Reds','YlOrBr','YlOrRd','hot']
	h_start = []
	h_end = []
	
	if highlights:
		h_start,h_end,_,_,_ = read_bedGraph(highFile,resolution,chromosome)
	
	rlength = len(files)-1 if compare and not pair else len(files)
	
	for exp in range(0,rlength):
		rowcounter=0
		
		if not tripleColumn:
			matrix,nums,tricks=read_HiCdata(files[exp],fileHeader,fileFooter,cleanNANs,smoothNoise,window,tadRange,plotInsulation,plotTadDomains,randomBins)
			end=len(matrix) if end == 0 else end
			if end > len(matrix): end=len(matrix)
			if start > len(matrix): start = 0
			size=end-start
			if exp == 0 : mlength = len(matrix)
			elif len(matrix) != mlength and not randomBins:
				print len(matrix), mlength
				print >>sys.stderr, 'unbalanced matrix size of '+files[exp]+' compared to '+files[0]+' ! matrix sizes should be equal'
				raise SystemExit
			
			matrix=matrix[start:end,start:end]
		else:
			if bedFile == '':
				print >>sys.stderr, 'an annotation bed file is required for triple-column sparse input.'
				raise SystemExit
			matrix,nums,tricks,clength=read_sparseHiCdata(files[exp],chromosome,bedFile,start,end,wholeGenome,smoothNoise,window,tadRange,plotInsulation,plotTadDomains,randomBins)
			if end > clength: end=clength
			if start > clength: start = 0
			end=clength if end == 0 else end
			size=end-start
	
		length = len(matrix)
		name=names[exp]	
		schr=chromosome.replace("chr","")
		
		
		''' MAIN matrix plotting '''
	
		ax1 = plt.subplot2grid((numOfrows, 4*len(files)), (0, exp*4), rowspan=4,colspan=4)

		ax1.set_title(('%s') % (name))
		if exp==0: 
			if not randomBins and not wholeGenome: ax1.set_ylabel('log2(interaction matrix) - %s Mb (resolution: %sKb)' % (chromosome , resolution/1000))
			elif randomBins: ax1.set_ylabel('log2(interaction matrix) - %s (Genomic Bins)' % (chromosome))
			elif wholeGenome: ax1.set_ylabel('')
			cmatrix = log2(pow(2, ceil(log2(max(matrix))/log2(2))))
			if matrixMax !=0: cmatrix = matrixMax
			if compare and not pair: matrix1=matrix
		
		if exp==1 and compare and not pair: matrix2=matrix
		if compare and pair: marray.append(matrix)
		 
		ax1.set_ylim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
		ax1.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			
		if not wholeGenome:
			if heatmapColor < 5:
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),cmap=plt.get_cmap(cmaps[heatmapColor]),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			elif heatmapColor == 3:
				cmap = plt.get_cmap(cmaps[heatmapColor])
				cmap.set_over('black')
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),cmap=cmap,origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			else:
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
		else:
			if heatmapColor < 5:
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),cmap=plt.get_cmap(cmaps[heatmapColor]),interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			elif heatmapColor == 3:
				cmap = plt.get_cmap(cmaps[heatmapColor])
				cmap.set_over('black')
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),cmap=cmap,origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			else:
				with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			plt.setp(ax1.get_xticklabels(), visible=False)
			
		if len(peakFiles) > 0:
			origin_x,origin_y,radius,colors = read_peakFile(peakFiles[exp],resolution,chromosome)
			for citem in range(0,len(origin_x)):
				if len(colors)==0: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor='black', linewidth=1, alpha=0.85)
				else: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor=colors[citem], linewidth=3, alpha=0.85)
				ax1.add_patch(circle)
				
		divider = make_axes_locatable(ax1)
		img.set_clim([0,cmatrix])
		
		if wholeGenome : plt.setp(ax1.get_yticklabels(), visible=False)
		ax1.get_yaxis().set_label_coords(-0.125,0.5) 
		
		if plotTadDomains and plotDomainTicks:
			ax1.set_xticks(tricks, minor=True)
			ax1.xaxis.grid(True,which='minor',linewidth=2)
		
		if h_start > 0:
			for item in range(0,len(h_start)):
				ax1.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
		
		rowcounter+=4
		ax1.get_xaxis().set_label_coords(0.5,-0.125)
		if numOfrows <= rowcounter and not randomBins and not wholeGenome: 
			cax = divider.append_axes("bottom", size="2.5%", pad=0.9)
			cbar = plt.colorbar(img, cax=cax, ticks=MultipleLocator(2.0), format="%.1f",orientation='horizontal',extendfrac='auto',spacing='uniform')
			plt.setp(ax1.get_xticklabels(), visible=True)
			ax1.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
		elif numOfrows <= rowcounter and randomBins: ax1.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
		elif numOfrows <= rowcounter and wholeGenome: ax1.set_xlabel('')
		else:
			cax = divider.append_axes("bottom", size="2.5%", pad=0.1)
			cbar = plt.colorbar(img, cax=cax, ticks=MultipleLocator(2.0), format="%.1f",orientation='horizontal',extendfrac='auto',spacing='uniform')
			plt.setp(ax1.get_xticklabels(), visible=False)
		
		''' Whole Genome matrix plotting '''
		
		if wholeGenome and numOfrows > rowcounter:
			print >>sys.stderr, 'Whole genome can be plotted only as matrix - this feature will be improved in future releases'
			raise SystemExit
		
		''' Triangular (Rotated Matrix) plotting '''
		
		if plotTriangular: 
		
			ax2 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=2,colspan=4,sharex=ax1)
			dst=ndimage.rotate(matrix,45,order=0,reshape=True,prefilter=False,cval=0)
			matrix=[];
			if not triangularHeight: height=length/5
			else: 
				if triangularHeight <= length: height = triangularHeight
				else: height = length/2
			
			ax2.set_ylim(start+length/2,start+length/2+height)
			ax2.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			ax2.set(adjustable='box-forced')
			if heatmapColor < 5:
				with np.errstate(divide='ignore'): img=ax2.imshow(log2(dst),origin="lower",cmap=plt.get_cmap(cmaps[heatmapColor]),interpolation="nearest",extent=(int(start or 1) - 0.5,\
															  	  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			else:
				with np.errstate(divide='ignore'): img=ax2.imshow(log2(dst),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
															  	  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			dst=[];
			img.set_clim([0,cmatrix-1])
			plt.setp(ax2.get_yticklabels(), visible=False)
			if exp==0: ax2.set_ylabel('Triangular')
			ax2.get_yaxis().set_label_coords(-0.125,0.5)
			if plotTadDomains and plotDomainTicks:
				ax2.set_xticks(tricks, minor=True)
				ax2.xaxis.grid(True,which='minor',linewidth=2)
			rowcounter+=2
			if numOfrows <= rowcounter and not randomBins: ax2.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax2.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			
			if h_start > 0:
				for item in range(0,len(h_start)):
					ax2.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
			
			ax2.spines['right'].set_visible(False)
			ax2.spines['top'].set_visible(False)
			ax2.spines['left'].set_visible(False)
			ax2.xaxis.set_ticks_position('bottom')
			plt.gca().yaxis.set_major_locator(plt.NullLocator())
			
		''' Random Bins matrix/triangular plotting '''
		
		if randomBins and numOfrows > rowcounter:
			print >>sys.stderr, 'Random bins data can be plotted only as matrix and triangular - this feature will be improved in future releases'
			raise SystemExit
		
		''' Gene plotting'''
		
		if len(plotGenes) > 0:
		
			ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=2,colspan=4,sharex=ax1)
			if exp==0: ax3.set_ylabel('Genes')
			ax3.get_yaxis().set_label_coords(-0.125,0.5)
			genes,trackCount,nearest = read_genes(plotGenes[0],resolution,chromosome,start,end)
			plength = (end-start)*float(resolution)/1000000
			
			if dark : gcolor='white';icolor='white'
			else: gcolor = '#3C3C8C';icolor='#0C0C78'
						
			for item in genes.keys():
								
				if len(genes[item])>2: #plot with introns
					gstart = float(item.split('-')[0])/resolution
					gend = float(item.split('-')[1])/resolution
					gtrack = genes[item][0]
					gestart = genes[item][3].split(',')
					geend = genes[item][4].split(',')
					
					ax3.plot([gstart,gend],[trackCount-gtrack+0.125,trackCount-gtrack+0.125],color=icolor, linewidth=0.5, zorder = -1)
					
					arrow = 5
					if genes[item][2]=='-': arrow=4
					
					if plength <= 30:
					
						for exon in range(0,len(geend)-1):
						
							if len(genes[item])>5:
								grect = Rectangle((float(gestart[exon])/resolution,trackCount-gtrack), (float(geend[exon])/resolution-float(gestart[exon])/resolution), 0.25, color=genes[item][5])
							else:
								grect = Rectangle((float(gestart[exon])/resolution,trackCount-gtrack), (float(geend[exon])/resolution-float(gestart[exon])/resolution), 0.25, color=gcolor)
						
						
							ax3.add_patch(grect)
							if exon < len(geend)-2:
								if (float(gestart[exon+1])/resolution-float(geend[exon])/resolution) > 0.5:
									if genes[item][2]=='-': ax3.plot(float(gestart[exon])/resolution+(float(gestart[exon+1])/resolution-float(geend[exon])/resolution)/2-0.125,trackCount-gtrack+0.125, marker=arrow, color=icolor, markersize=1.25)
									else: ax3.plot(float(gestart[exon])/resolution+(float(gestart[exon+1])/resolution-float(geend[exon])/resolution)/2+0.125,trackCount-gtrack+0.125, marker=arrow, color=icolor, markersize=1.25)					
												
					else:
							
						if len(genes[item])>5: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=genes[item][5])
						else: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=gcolor)
						ax3.add_patch(rect)			
								
				else: # simple plotting
										
					gstart = float(item.split('-')[0])/resolution
					gend = float(item.split('-')[1])/resolution
					gtrack = genes[item][0]
					
					if len(genes[item])>5: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=genes[item][5])
					else: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=gcolor)
					ax3.add_patch(rect)
									
 				
 				if plength <= 30 and geneLabels: # also consider the gene density
					
					#optimize the font size	
 				
 					gindex = nearest[gtrack].index(int(item.split('-')[1]))
 					upgene = nearest[gtrack][gindex-1]
 					if gindex < len(nearest[gtrack])-1: downgene = nearest[gtrack][gindex+1] 
 					else: downgene = upgene
 
					if plength <= 2 or plength < 1: plength=1
					elif plength <= 4: plength = 2
					else : plength/1.5
					gdist = min(abs(nearest[gtrack][gindex]-upgene),abs(nearest[gtrack][gindex]-downgene))
					if len(nearest[gtrack])==1: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=4.5/plength)
					elif float(gdist)/resolution >= 2 and len(genes[item][1])<=6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=4.5/plength)
					elif float(gdist)/resolution >= 2 and len(genes[item][1])>6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=3/plength)
					elif float(gdist)/resolution < 2 and float(gdist)/resolution > 1 and len(genes[item][1])<=6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=2.5/plength)
					elif float(gdist)/resolution < 2 and float(gdist)/resolution >1 and len(genes[item][1])>6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=2/plength)
					elif float(gdist)/resolution <= 1 and float(gdist)/resolution >= 0.25: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=1.8)
					else: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=1)
				
					#if lblstndrd == 1:
 						#if len(genes[item][1])>5: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=0.5)
 						#else : ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=3)	
			
			ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			ax3.set_ylim(0,trackCount+1)
			plt.setp(ax3.get_yticklabels(), visible=False)
			if plotTadDomains and plotDomainTicks:
				ax3.set_xticks(tricks, minor=True)
				ax3.xaxis.grid(True,which='minor')
			
			if h_start > 0:
				for item in range(0,len(h_start)):
					ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
			
			ax3.spines['right'].set_visible(False)
			ax3.spines['left'].set_visible(False)
			ax3.spines['top'].set_visible(False)
			ax3.tick_params(left="off")
			ax3.tick_params(right="off")
			ax3.xaxis.set_ticks_position('bottom')
			
			rowcounter+=2
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			
		
		''' Histogram plotting '''
		
		if len(histograms)>0:
			if not superImpose: # improve this part
				for x in range(0,len(histograms[0].split(','))):
					ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
					ax3.get_yaxis().set_label_coords(-0.125,0.5)
					x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(histograms[exp].split(',')[x],resolution,chromosome)
					
					if dark: ax3.plot(x_comps,y_comps,color='white')
					else: ax3.plot(x_comps,y_comps,color='black')
					
					if exp==0: 
						ystart,yend = where(start,end,x_comps)
						ymin = min(y_comps[ystart:yend])+ min(y_comps[ystart:yend])/10 if min(y_comps[ystart:yend]) < 0 else min(y_comps[ystart:yend])-min(y_comps[ystart:yend])/10
						yminlims.append(ymin)
						if len(histMax)==0:
							ax3.set_ylim(ymin,max(y_comps[ystart:yend])+max(y_comps[ystart:yend])/10)
							ymaxlims.append(max(y_comps[ystart:yend]))
							#print ymin, max(y_comps[ystart:yend])+max(y_comps[ystart:yend])/10
						else:
							ax3.set_ylim(ymin,int(histMax[exp].split(',')[x])+int(histMax[exp].split(',')[x])/10)
						ax3.set_ylabel(histLabels[exp].split(',')[x])
					else:
						if len(histMax)==0:
							ax3.set_ylim(yminlims[x],ymaxlims[x]+ymaxlims[x]/10)
						else:
							ax3.set_ylim(ymin,int(histMax[exp].split(',')[x])+int(histMax[exp].split(',')[x])/10)
			
					ax3.locator_params(axis='y',tight=False, nbins=3)
					ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			
			
					if len(fillHist) > 0 and int(fillHist[exp].split(',')[x])==1:
						comps2=array(y_comps)
						if len(histColors)>0:
							if histColors[exp].split(',')[x] != '':
								ax3.fill_between(x_comps, comps2,0, color='#'+histColors[exp].split(',')[x], interpolate=True)
								if ymin < 0: 
									with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='#'+histColors[exp].split(',')[x], interpolate=True)
									with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
							else:
								ax3.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
								if ymin < 0: 
									with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
									with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
						else:
							ax3.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
							if ymin < 0: 
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
						
					x_comps=[];x_comps2=[];y_comps=[];colors==[];	
					if plotTadDomains and plotDomainTicks:
						ax3.set_xticks(tricks, minor=True)
						ax3.xaxis.grid(True,which='minor')
			
					if h_start > 0:
						for item in range(0,len(h_start)):
							ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
			
					if spine > 0:
						ax3.spines['right'].set_visible(False)
						ax3.spines['top'].set_visible(False)
						ax3.xaxis.set_ticks_position('bottom')
						ax3.yaxis.set_ticks_position('left')
			
					rowcounter+=1
				
			else:
				
				hfirst = histograms[exp].split(',')[0]
				hsecond = histograms[exp].split(',')[1]
				
				ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
				ax3.get_yaxis().set_label_coords(-0.125,0.5)
					
				x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(hfirst,resolution,chromosome)
				
				if dark : firstcolor='white'
				else: firstcolor = 'black'
				
				if len(histColors)>0:
					if histColors[exp].split(',')[0] != '': ax3.plot(x_comps,y_comps,color='#'+histColors[exp].split(',')[0],alpha=0.5);firstcolor='#'+histColors[exp].split(',')[0]
				elif dark: x3.plot(x_comps,y_comps,color='white',alpha=0.5)
				else: ax3.plot(x_comps,y_comps,color='black',alpha=0.5)
				
				ax3.tick_params(axis='y', colors=firstcolor)
				
				if exp==0: 
					ystart,yend = where(start,end,x_comps)
					ymin = min(y_comps[ystart:yend])+ min(y_comps[ystart:yend])/10 if min(y_comps[ystart:yend]) < 0 else min(y_comps[ystart:yend])-min(y_comps[ystart:yend])/10
					yminlims.append(ymin)
					if len(histMax)==0:
						ax3.set_ylim(ymin,max(y_comps[ystart:yend])+max(y_comps[ystart:yend])/10)
						ymaxlims.append(max(y_comps[ystart:yend]))
					else:
						ax3.set_ylim(ymin,int(histMax[exp].split(',')[0])+int(histMax[exp].split(',')[0])/10)
					ax3.set_ylabel(histLabels[exp].split(',')[0],color=firstcolor)
				else:
					if len(histMax)==0:
						ax3.set_ylim(yminlims[0],ymaxlims[0]+ymaxlims[0]/10)
					else:
						ax3.set_ylim(ymin,int(histMax[exp].split(',')[0])+int(histMax[exp].split(',')[0])/10)
				
				
				if len(fillHist) > 0 and int(fillHist[exp].split(',')[0])==1:
					comps2=array(y_comps)
					if len(histColors)>0:
						if histColors[exp].split(',')[0] != '':
							ax3.fill_between(x_comps, comps2,0, color='#'+histColors[exp].split(',')[0], interpolate=True)
							if ymin < 0: 
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='#'+histColors[exp].split(',')[0], interpolate=True)
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
						else:
							ax3.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
							if ymin < 0: 
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
								with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
					else:
						ax3.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
						if ymin < 0: 
							with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
							with np.errstate(all='ignore'):ax3.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)

				
				ax3.locator_params(axis='y',tight=False, nbins=3)
				ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
					
				ax3a = ax3.twinx() # second plot
				subcolor = 'red'
				x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(hsecond,resolution,chromosome)
				if len(histColors)>0:
					if histColors[exp].split(',')[1] != '': ax3a.plot(x_comps,y_comps,color='#'+histColors[exp].split(',')[1]); subcolor= '#'+histColors[exp].split(',')[1]
				else: ax3a.plot(x_comps,y_comps,color='red')
				
				ax3a.patch.set_visible(False)
				if exp==0: 
					ystart,yend = where(start,end,x_comps)
					ymin = min(y_comps[ystart:yend])+ min(y_comps[ystart:yend])/10 if min(y_comps[ystart:yend]) < 0 else min(y_comps[ystart:yend])-min(y_comps[ystart:yend])/10
					yminlims.append(ymin)
					if len(histMax)==0:
						ax3a.set_ylim(ymin,max(y_comps[ystart:yend])+max(y_comps[ystart:yend])/10)
						ymaxlims.append(max(y_comps[ystart:yend]))
					else:
						ax3a.set_ylim(ymin,int(histMax[exp].split(',')[1])+int(histMax[exp].split(',')[1])/10)
					ax3a.text(int(start or 1)+length-length/6,max(y_comps[ystart:yend])-max(y_comps[ystart:yend])/10,histLabels[exp].split(',')[1],color=subcolor)
					#ax3a.set_ylabel(histLabels[exp].split(',')[1])
				else:
					if len(histMax)==0:
						ax3a.set_ylim(yminlims[1],ymaxlims[1]+ymaxlims[1]/10)
					else:
						ax3a.set_ylim(ymin,int(histMax[exp].split(',')[1])+int(histMax[exp].split(',')[1])/10)
				
				if len(fillHist) > 0 and int(fillHist[exp].split(',')[1])==1:
					comps2=array(y_comps)
					if len(histColors)>0:
						if histColors[exp].split(',')[1] != '':
							ax3a.fill_between(x_comps, comps2,0, color='#'+histColors[exp].split(',')[1], interpolate=True)
							if ymin < 0: 
								with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2,0, comps2>0, color='#'+histColors[exp].split(',')[1], interpolate=True)
								with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
						else:
							ax3.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
							if ymin < 0: 
								with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
								with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
					else:
						ax3a.fill_between(x_comps, comps2,0, color='gray', interpolate=True)
						if ymin < 0: 
							with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2,0, comps2>0, color='gray', interpolate=True)
							with np.errstate(all='ignore'):ax3a.fill_between(x_comps, comps2, 0, where=comps2<0, color='black', interpolate=True)
				
				ax3a.locator_params(axis='y',tight=False, nbins=3)
				ax3a.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
				ax3a.tick_params(axis='y', colors=subcolor)
				x_comps=[];x_comps2=[];y_comps=[];colors==[];	
				
				if plotTadDomains and plotDomainTicks:
					ax3.set_xticks(tricks, minor=True)
					ax3.xaxis.grid(True,which='minor')
		
				if h_start > 0:
					for item in range(0,len(h_start)):
						ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
		
				if spine > 0:
					ax3.spines['right'].set_visible(False)
					ax3.spines['top'].set_visible(False)
					ax3.xaxis.set_ticks_position('bottom')
					ax3.yaxis.set_ticks_position('left')
					ax3a.spines['left'].set_visible(False)
					ax3a.spines['top'].set_visible(False)
					ax3a.xaxis.set_ticks_position('bottom')
					ax3a.yaxis.set_ticks_position('right')

				rowcounter+=1
				
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
		
		
		''' Bar plotting '''
		
		if len(barPlots)>0: 
			for x in range(0,len(barPlots[0].split(','))):
	
				ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
				if exp==0: ax3.set_ylabel(barLabels[exp].split(',')[x])
				ax3.get_yaxis().set_label_coords(-0.125,0.5)
				x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(barPlots[exp].split(',')[x],resolution,chromosome)
				
				if len(barMax)==0: hMax = max(y_comps)
				else: hMax = float(barMax[exp].split(',')[x]) #need to implement length check
				
				for item in range(0,len(x_comps)):
					#if x_comps[item]>=start and x_comps[item]<=end:
					if len(barMax)>0 and y_comps[item]>float(barMax[exp].split(',')[x]): y_comps[item]=float(barMax[exp].split(',')[x])
					if len(barColors)==0 and len(colors)==0: rect = Rectangle((x_comps[item],0.0), (x_comps2[item]-x_comps[item]), y_comps[item], color='#0099FF',alpha=y_comps[item]/hMax)
					elif len(colors)>0: rect = Rectangle((x_comps[item],0.0), (x_comps2[item]-x_comps[item]),  y_comps[item], color=colors[item])
					elif len(barColors)>0: rect = Rectangle((x_comps[item],0.0), (x_comps2[item]-x_comps[item]),  y_comps[item], color='#'+barColors[exp].split(',')[x],alpha=y_comps[item]/hMax)
					ax3.add_patch(rect)
				x_comps=[];x_comps2=[];y_comps=[];colors==[];
				ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
				ax3.set_ylim(0,hMax+hMax/10)
				ax3.locator_params(axis='y',tight=False, nbins=4)
				if plotTadDomains and plotDomainTicks:
					ax3.set_xticks(tricks, minor=True)
					ax3.xaxis.grid(True,which='minor')
				
				if h_start > 0:
					for item in range(0,len(h_start)):
						ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
				
				if spine > 0:
					ax3.spines['right'].set_visible(False)
					ax3.spines['top'].set_visible(False)
					ax3.xaxis.set_ticks_position('bottom')
					ax3.yaxis.set_ticks_position('left')
				
				rowcounter+=1
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))

		
				
		''' Tile plotting '''
		
		if len(tilePlots)>0: 
			for x in range(0,len(tilePlots[0].split(','))):
	
				ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
				if exp==0: ax3.set_ylabel(tileLabels[exp].split(',')[x])
				ax3.get_yaxis().set_label_coords(-0.125,0.5)
				x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(tilePlots[exp].split(',')[x],resolution,chromosome)
				for item in range(0,len(x_comps)):
					if len(tileColors)==0 and len(colors)==0: rect = Rectangle((x_comps[item],0.35), (x_comps2[item]-x_comps[item]), 0.25, color='#0099FF')
					elif len(colors)>0: rect = Rectangle((x_comps[item],0.35), (x_comps2[item]-x_comps[item]), 0.25, color=colors[item])
					elif len(tileColors)>0: rect = Rectangle((x_comps[item],0.35), (x_comps2[item]-x_comps[item]), 0.25, color='#'+tileColors[exp].split(',')[x])
					if len(texts) > 0 and tileText: ax3.text(x_comps[item]-1, 0.75, texts[item], fontsize=10)
					ax3.add_patch(rect)
				x_comps=[];x_comps2=[];y_comps=[];colors==[];
				ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
				ax3.set_ylim(0,1)
				plt.setp(ax3.get_yticklabels(), visible=False)
				if plotTadDomains and plotDomainTicks:
					ax3.set_xticks(tricks, minor=True)
					ax3.xaxis.grid(True,which='minor')
				
				if h_start > 0:
					for item in range(0,len(h_start)):
						ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
				
				if spine > 0:
					ax3.spines['right'].set_visible(False)
					ax3.spines['top'].set_visible(False)
					ax3.xaxis.set_ticks_position('bottom')
					ax3.yaxis.set_ticks_position('left')
				
				rowcounter+=1
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
		
		
		''' Arc plotting '''
		
		if len(arcPlots)>0: 
			for x in range(0,len(arcPlots[0].split(','))):

				ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
				if exp==0: ax3.set_ylabel(arcLabels[exp].split(',')[x])
				ax3.get_yaxis().set_label_coords(-0.125,0.5)
				x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(arcPlots[exp].split(',')[x],resolution,chromosome)
				ymax = 0
				
				for item in range(0,len(x_comps)):
				
					center = x_comps[item]+(x_comps2[item]-x_comps[item])/2
					rad = (x_comps2[item]-x_comps[item])/2
					pts = get_ellipse_coords(a=rad, b=1.0, x=center, k=1./8)
					if dark and len(arcColors)==0 and len(colors)==0: ax3.plot(pts[:,0], pts[:,1],c='#cecece')
					elif len(arcColors)==0 and len(colors)==0: ax3.plot(pts[:,0], pts[:,1],c='black')
					elif len(colors)>0: ax3.fill_between(pts[:,0], pts[:,1],0, color=colors[item], interpolate=True, alpha=0.35)
					elif len(arcColors)>0: ax3.fill_between(pts[:,0], pts[:,1],0, color='#'+arcColors[exp].split(',')[x], interpolate=True, alpha=0.35)
				
				x_comps=[];x_comps2=[];y_comps=[];colors==[];
				ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
				ax3.set_ylim(0,1)
				plt.setp(ax3.get_yticklabels(), visible=False)
				if plotTadDomains and plotDomainTicks:
					ax3.set_xticks(tricks, minor=True)
					ax3.xaxis.grid(True,which='minor')
				
				if h_start > 0:
					for item in range(0,len(h_start)):
						ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
				
				if spine > 0:
					ax3.spines['right'].set_visible(False)
					ax3.spines['top'].set_visible(False)
					ax3.xaxis.set_ticks_position('bottom')
					ax3.yaxis.set_ticks_position('left')
				
				rowcounter+=1
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
		
		''' Epilogos plotting '''
		
		if len(epiLogos)>0:
			if imputed:
				obs_colors={'1':'#ff0000','2':'#ff4500','3':'#ff4500','4':'#ff4500','5':'#008000','6':'#008000','7':'#008000','8':'#009600','9':'#c2e105','10':'#c2e105','11':'#c2e105','12':'#c2e105','13':'#ffc34d','14':'#ffc34d','15':'#ffc34d','16':'#ffff00','17':'#ffff00','18':'#ffff00','19':'#ffff66','20':'#66cdaa','21':'#8a91d0','22':'#e6b8b7','23':'#7030a0','24':'#808080','25':'#ffffff'} 
			else:
				obs_colors={'1':'#ff0000','2':'#ff4500','3':'#32cd32','4':'#008000','5':'#006400','6':'#c2e105','7':'#ffff00','8':'#66cdaa','9':'#8a91d0','10':'#cd5c5c','11':'#e9967a','12':'#bdb76b','13':'#808080','14':'#c0c0c0','15':'#ffffff'}
			
			ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
			if exp==0: ax3.set_ylabel('Epilogos')
			ax3.get_yaxis().set_label_coords(-0.125,0.5)
			x_comps,y_dict = read_epilogos(epiLogos,resolution,chromosome,start,end)				
			for state in y_dict.keys():
				ax3.plot(x_comps,y_dict[state],color=obs_colors[state],alpha=0.75,linewidth = 0.8)
				if imputed:
					if state not in ['9','10','11','12','21']:
						ax3.plot(x_comps,y_dict[state],color=obs_colors[state],alpha=0.45,linewidth = 0.8)
					elif state =='21':
						ax3.plot(x_comps,y_dict[state],color=obs_colors[state],alpha=0.75,linewidth = 0.8)
			ax3.set_ylim(1,7)
			ax3.locator_params(axis='y',tight=False, nbins=3)
			ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			ax3.set_axis_bgcolor('black')
			
			if spine > 0:
				ax3.spines['right'].set_visible(False)
				ax3.spines['top'].set_visible(False)
				ax3.xaxis.set_ticks_position('bottom')
				ax3.yaxis.set_ticks_position('left')
			
			rowcounter+=1
			if numOfrows <= rowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)

		
		
		''' Insulation Scores '''
		
		if plotInsulation: 
		
			ax4 = plt.subplot2grid((numOfrows,4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
			if exp==0 and tripleColumn:
				ins_score = max(nums)+max(nums)/5
			elif exp==0: 
				ax4.set_ylabel('Insulation')
				ins_score = max(nums[start:end])+max(nums[start:end])/5
					
			ax4.get_yaxis().set_label_coords(-0.125,0.5)
			ax4.locator_params(axis='y',tight=False, nbins=3)		
			ax4.set_ylim(0,ins_score)
			ax4.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			if not tripleColumn:
				if dark: ax4.plot(range(0,len(nums)),nums,'white')
				else: ax4.plot(range(0,len(nums)),nums,'black')
				ax4.fill_between(range(0,len(nums)),nums,0,color='0.8')
			else:
				if dark: ax4.plot(np.arange(start,end),nums,'white')
				else: ax4.plot(np.arange(start,end),nums,'black')
				ax4.fill_between(np.arange(start,end),nums,0,color='0.8')
				
			if plotTadDomains and plotDomainTicks:
				ax4.set_xticks(tricks, minor=True)
				ax4.xaxis.grid(True,which='minor')
			
			if h_start > 0:
				for item in range(0,len(h_start)):
					ax4.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
				
			if spine > 0:
				ax4.spines['right'].set_visible(False)
				ax4.spines['top'].set_visible(False)
				ax4.xaxis.set_ticks_position('bottom')
				ax4.yaxis.set_ticks_position('left')
					
			rowcounter+=1
			if numOfrows <= rowcounter and not randomBins: ax4.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000)) 
			elif numOfrows <= rowcounter and randomBins: ax4.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			
		'''TAD plotings - determined by insulation score'''
		
		if plotTadDomains: 
	
			ax5 = plt.subplot2grid((numOfrows,4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
			
			for item in range(0,len(tricks)-1):
				
				tcolor='darkkhaki'
				if len(domColors)>0: tcolor='#'+domColors[exp]
						
				if plotDomainsAsBars:
					if not plotPublishedTadDomains: p = Rectangle((tricks[item],0.2), (tricks[item+1]-tricks[item]), 0.25, color=tcolor,alpha=0.75)
					else: p = Rectangle((tricks[item],0.1), (tricks[item+1]-tricks[item]), 0.15, color=tcolor,alpha=0.75)
				else:
					if not tripleColumn:
						pts= np.array([[tricks[item],0],[tricks[item+1],0],[floor((tricks[item]+tricks[item+1])/2),0.75]])
						p = Polygon(pts, closed=True,color=tcolor,alpha=max(nums[tricks[item]:tricks[item+1]])/max(nums))
					else:
						pts= np.array([[tricks[item],0],[tricks[item+1],0],[floor((tricks[item]+tricks[item+1])/2),0.75]])
						p = Polygon(pts, closed=True,color=tcolor)
				
				if not tripleColumn and sum(nums[slice(tricks[item],tricks[item+1])]) > np.percentile(np.array(nums),75):
					ax5.add_patch(p)
				elif tripleColumn and sum(nums[slice(tricks[item]-start,tricks[item+1]-start)]) > np.percentile(np.array(nums),75):
					ax5.add_patch(p)
					
			if plotPublishedTadDomains:
				## adding TAD domain predictions from Dixon et al. Nature 2009
				if publishedTadDomainOrganism:
					fone=open('data/IMR90_domains_hg19.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.3), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
								p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
							ax5.add_patch(p)
					fone=open('data/hESC_domains_hg19.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.5), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
								p = Polygon(pts, closed=True,color='steelblue',alpha=0.5)
							ax5.add_patch(p)
					ax5.set_title('Khaki:%s - Blue:hES - Red:IMR90' % (name),fontsize=8)
				else:
					fone=open('data/mCortex_domains_mm9.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.3), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
								p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
							ax5.add_patch(p)
					fone=open('data/mES_domains_mm9.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.5), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
								p = Polygon(pts, closed=True,color='steelblue',alpha=0.5)
							ax5.add_patch(p)
					ax5.set_title('Khaki:%s - Blue:mES - Red:Cortex' % (name),fontsize=8)
			else:
				ax5.set_title('Khaki:%s' % (name))
				
			ax5.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			if exp==0: ax5.set_ylabel("Domains")
			ax5.locator_params(axis='y',tight=False, nbins=3)
			ax5.get_yaxis().set_label_coords(-0.125,0.5)
			plt.setp(ax5.get_yticklabels(), visible=False)
			if numOfrows <= rowcounter and not plotCustomDomains and not randomBins: ax5.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and not plotCustomDomains and randomBins: ax5.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			ax5.set_ylim(0,0.75)
			if spine > 0:
				ax5.spines['right'].set_visible(False)
				ax5.spines['top'].set_visible(False)
				ax5.xaxis.set_ticks_position('bottom')
				ax5.yaxis.set_ticks_position('left')
			rowcounter+=1
		
		'''TAD plotings - custom domains'''
		
		if plotCustomDomains and not plotTadDomains:
			x_comps,x_comps2,y_comps,colors,texts = read_bedGraph(customDomainsFile[exp],resolution,chromosome)
			ax5 = plt.subplot2grid((numOfrows,4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
			
			for item in range(0,len(x_comps)):
				tcolor='darkkhaki'
				if len(colors)>0: tcolor=colors[item]
				elif len(domColors)>0: tcolor='#'+domColors[exp]
				
				if plotDomainsAsBars:
					if not plotPublishedTadDomains: p = Rectangle((x_comps[item],0.2), (x_comps2[item]-x_comps[item]), 0.25, color=tcolor,alpha=0.75)
					else: p = Rectangle((x_comps[item],0.1), (x_comps2[item]-x_comps[item]), 0.15, color=tcolor,alpha=0.75)
				else:
					pts= np.array([[x_comps[item],0],[x_comps2[item],0],[floor((x_comps[item]+x_comps2[item])/2),0.75]])
					p = Polygon(pts, closed=True,color=tcolor,alpha=0.85)
				ax5.add_patch(p)
			
			if plotPublishedTadDomains:
				## adding TAD domain predictions from Dixon et al. Nature 2009 - genome assemblies hg19,mm9
				if publishedTadDomainOrganism:
					fone=open('data/IMR90_domains_hg19.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.3), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
								p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
							ax5.add_patch(p)
					fone=open('data/hESC_domains_hg19.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.5), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
								p = Polygon(pts, closed=True,color='steelblue',alpha=0.4)
							ax5.add_patch(p)
					ax5.set_title('Khaki:%s - Blue:hES - Red:IMR90' % (name),fontsize=8)
				else:
					fone=open('data/mES_domains_mm9.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.3), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
								p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
							ax5.add_patch(p)
					fone=open('data/mCortex_domains_mm9.bed','r')
					for line in fone.xreadlines():
						tags = line.strip().split("\t")
						if tags[0]==chromosome:
							Tstart = int(tags[1])/resolution
							Tend = int(tags[2])/resolution
							if plotDomainsAsBars:
								p = Rectangle((Tstart,0.5), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
							else:
								pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
								p = Polygon(pts, closed=True,color='steelblue',alpha=0.4)
							ax5.add_patch(p)
					ax5.set_title('Khaki:%s - Red:mES - Blue:Cortex' % (name),fontsize=8)
			else:
				ax5.set_title('Khaki:%s' % (name),fontsize=8)
				
			ax5.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			if exp==0: ax5.set_ylabel("Domains")
			ax5.locator_params(axis='y',tight=False, nbins=3)
			ax5.get_yaxis().set_label_coords(-0.125,0.5)
			plt.setp(ax5.get_yticklabels(), visible=False)
			if numOfrows <= rowcounter and not randomBins : ax5.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins: ax5.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			ax5.set_ylim(0,0.75)
			if spine > 0:
				ax5.spines['right'].set_visible(False)
				ax5.spines['top'].set_visible(False)
				ax5.xaxis.set_ticks_position('bottom')
				ax5.yaxis.set_ticks_position('left')
			rowcounter+=1
		elif plotPublishedTadDomains and not plotTadDomains:
			ax5 = plt.subplot2grid((numOfrows,4*len(files)), (rowcounter, exp*4), rowspan=1,colspan=4,sharex=ax1)
			if publishedTadDomainOrganism:
				fone=open('data/IMR90_domains_hg19.bed','r')
				for line in fone.xreadlines():
					tags = line.strip().split("\t")
					if tags[0]==chromosome:
						Tstart = int(tags[1])/resolution
						Tend = int(tags[2])/resolution
						if plotDomainsAsBars:
							p = Rectangle((Tstart,0.2), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
						else:
							pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
							p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
						ax5.add_patch(p)
				fone=open('data/hESC_domains_hg19.bed','r')
				for line in fone.xreadlines():
					tags = line.strip().split("\t")
					if tags[0]==chromosome:
						Tstart = int(tags[1])/resolution
						Tend = int(tags[2])/resolution
						if plotDomainsAsBars:
							p = Rectangle((Tstart,0.4), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
						else:
							pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
							p = Polygon(pts, closed=True,color='steelblue',alpha=0.4)
						ax5.add_patch(p)
				ax5.set_title('Blue:hES - Red:IMR90',fontsize=8)
			else:
				fone=open('data/mCortex_domains_mm9.bed','r')
				for line in fone.xreadlines():
					tags = line.strip().split("\t")
					if tags[0]==chromosome:
						Tstart = int(tags[1])/resolution
						Tend = int(tags[2])/resolution
						if plotDomainsAsBars:
							p = Rectangle((Tstart,0.2), (Tend-Tstart), 0.15, color='salmon',alpha=0.75)
						else:
							pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.25]])
							p = Polygon(pts, closed=True,color='salmon',alpha=0.5)
						ax5.add_patch(p)
				fone=open('data/mES_domains_mm9.bed','r')
				for line in fone.xreadlines():
					tags = line.strip().split("\t")
					if tags[0]==chromosome:
						Tstart = int(tags[1])/resolution
						Tend = int(tags[2])/resolution
						if plotDomainsAsBars:
							p = Rectangle((Tstart,0.4), (Tend-Tstart), 0.15, color='steelblue',alpha=0.75)
						else:
							pts= np.array([[Tstart,0],[Tend,0],[floor((Tstart+Tend)/2),0.4]])
							p = Polygon(pts, closed=True,color='steelblue',alpha=0.4)
						ax5.add_patch(p)
				ax5.set_title('Blue:mES - Red:Cortex',fontsize=8)
				
			ax5.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			if exp==0: ax5.set_ylabel("Domains")
			ax5.locator_params(axis='y',tight=False, nbins=3)
			ax5.get_yaxis().set_label_coords(-0.125,0.5)
			plt.setp(ax5.get_yticklabels(), visible=False)
			if numOfrows <= rowcounter and not randomBins : ax5.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= rowcounter and randomBins : ax5.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			ax5.set_ylim(0,0.75)
			if spine > 0:
				ax5.spines['right'].set_visible(False)
				ax5.spines['top'].set_visible(False)
				ax5.xaxis.set_ticks_position('bottom')
				ax5.yaxis.set_ticks_position('left')
			rowcounter+=1

		if not randomBins:
			ticks= ax1.get_xticks().tolist()
			if resolution*(end-start)<=500000: 
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,3)
			elif resolution*(end-start)<=1500000:
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,2)
			else: 
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,1) 
			ax1.set_xticklabels(ticks)
			ax1.set_yticklabels(ticks)

	# Single comparisons
	if compare and not pair:
		crowcounter = 4
		if len(compareSm)>0 :
		
			slow = int(compareSm.split(',')[0])
			shigh = int(compareSm.split(',')[1])
			matrix1[(matrix1>=slow) & (matrix1<=shigh)]=1
			matrix2[(matrix2>=slow) & (matrix2<=shigh)]=1

		with np.errstate(divide='ignore',invalid='ignore'): matrix=matrix1/matrix2
		#matrix[np.logical_and(matrix>=0.5, matrix<=1)]=1
		ax1 = plt.subplot2grid((numOfrows, 4*len(files)), (0, (exp+1)*4), rowspan=4,colspan=4)
		with np.errstate(divide='ignore'): img=ax1.imshow(log2(matrix),cmap=plt.get_cmap("bwr"),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
		
		ax1.set_ylim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
		ax1.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
		ax1.set_title(('log2(%s / %s)') % (names[0],names[1]))
		
		if len(peakFiles) > 0:
			origin_x,origin_y,radius,colors = read_peakFile(peakFiles[exp],resolution,chromosome)
			for citem in range(0,len(origin_x)):
				if len(colors)==0: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor='black', linewidth=1, alpha=0.85)
				else: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor=colors[citem], linewidth=3, alpha=0.85)
				ax1.add_patch(circle)
				
		divider = make_axes_locatable(ax1)
		if len(compareEx)>0 : clow = int(compareEx.split(',')[0])*-1;chigh=int(compareEx.split(',')[1]);img.set_clim([clow,chigh])
		else: img.set_clim([-4,4])
				
		if wholeGenome : plt.setp(ax1.get_yticklabels(), visible=False)
		ax1.get_yaxis().set_label_coords(-0.125,0.5) 
		if plotTadDomains and plotDomainTicks:
			ax1.set_xticks(tricks, minor=True)
			ax1.xaxis.grid(True,which='minor',linewidth=2)
		
		if h_start > 0:
			for item in range(0,len(h_start)):
				ax1.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
		
		ax1.get_xaxis().set_label_coords(0.5,-0.125)
		if numOfrows == 4 and not randomBins:
			cax = divider.append_axes("bottom", size="2.5%", pad=0.9)
			cbar = plt.colorbar(img, cax=cax, ticks=MultipleLocator(2.0), format="%.1f",orientation='horizontal',extendfrac='auto',spacing='uniform')
			plt.setp(ax1.get_xticklabels(), visible=True)
			ax1.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))	
		else:
			cax = divider.append_axes("bottom", size="2.5%", pad=0.1)
			cbar = plt.colorbar(img, cax=cax, ticks=MultipleLocator(2.0), format="%.1f",orientation='horizontal',extendfrac='auto',spacing='uniform')
			plt.setp(ax1.get_xticklabels(), visible=False)
		
		if not randomBins:
			ticks= ax1.get_xticks().tolist()
			if resolution*(end-start)<=500000: 
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,3)
			elif resolution*(end-start)<=1500000:
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,2)
			else: 
				for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,1) 
			ax1.set_xticklabels(ticks)
			ax1.set_yticklabels(ticks)
		
		if plotTriangular: 
			
			ax2 = plt.subplot2grid((numOfrows, 4*len(files)), (crowcounter, (exp+1)*4), rowspan=2,colspan=4,sharex=ax1)
			dst=ndimage.rotate(matrix,45,order=0,reshape=True,prefilter=False,cval=0)
			matrix=[];
			if not triangularHeight: height=length/5
			else: 
				if triangularHeight <= length: height = triangularHeight
				else: height = length/2
			
			ax2.set_ylim(start+length/2,start+length/2+height)
			ax2.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			ax2.set(adjustable='box-forced')
			with np.errstate(divide='ignore'): img=ax2.imshow(log2(dst),cmap=plt.get_cmap("bwr"),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
														  		  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
			dst=[];
			if len(compareEx)>0 : clow = int(compareEx.split(',')[0])*-1;chigh=int(compareEx.split(',')[1]);img.set_clim([clow,chigh])
			else: img.set_clim([-4,4])

			plt.setp(ax2.get_yticklabels(), visible=False)
			if exp==0: ax2.set_ylabel('Triangular')
			ax2.get_yaxis().set_label_coords(-0.125,0.5)
			if plotTadDomains and plotDomainTicks:
				ax2.set_xticks(tricks, minor=True)
				ax2.xaxis.grid(True,which='minor',linewidth=2)
			crowcounter+=2
			if numOfrows <= crowcounter and not randomBins: ax2.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= crowcounter and randomBins: ax2.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
			
			if h_start > 0:
				for item in range(0,len(h_start)):
					ax2.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
			
			ax2.spines['right'].set_visible(False)
			ax2.spines['top'].set_visible(False)
			ax2.spines['left'].set_visible(False)
			ax2.xaxis.set_ticks_position('bottom')
			plt.gca().yaxis.set_major_locator(plt.NullLocator())
		
		if len(plotGenes) > 0:
		
			ax3 = plt.subplot2grid((numOfrows, 4*len(files)), (crowcounter, (exp+1)*4), rowspan=2,colspan=4,sharex=ax1)
			if exp==0: ax3.set_ylabel('Genes')
			ax3.get_yaxis().set_label_coords(-0.125,0.5)
			genes,trackCount,nearest = read_genes(plotGenes[0],resolution,chromosome,start,end)
			plength = (end-start)*float(resolution)/1000000
						
			for item in genes.keys():
								
				if len(genes[item])>2: #plot with introns
					gstart = float(item.split('-')[0])/resolution
					gend = float(item.split('-')[1])/resolution
					gtrack = genes[item][0]
					gestart = genes[item][3].split(',')
					geend = genes[item][4].split(',')
					
					ax3.plot([gstart,gend],[trackCount-gtrack+0.125,trackCount-gtrack+0.125],color=icolor, linewidth=0.5, zorder = -1)
					
					arrow = 5
					if genes[item][2]=='-': arrow=4
					
					if plength <= 30:
					
						for exon in range(0,len(geend)-1):
						
							if len(genes[item])>5:
								grect = Rectangle((float(gestart[exon])/resolution,trackCount-gtrack), (float(geend[exon])/resolution-float(gestart[exon])/resolution), 0.25, color=genes[item][5])
							else:
								grect = Rectangle((float(gestart[exon])/resolution,trackCount-gtrack), (float(geend[exon])/resolution-float(gestart[exon])/resolution), 0.25, color=gcolor)
						
						
							ax3.add_patch(grect)
							if exon < len(geend)-2:
								if (float(gestart[exon+1])/resolution-float(geend[exon])/resolution) > 0.5:
									if genes[item][2]=='-': ax3.plot(float(gestart[exon])/resolution+(float(gestart[exon+1])/resolution-float(geend[exon])/resolution)/2-0.125,trackCount-gtrack+0.125, marker=arrow, color=icolor, markersize=1.25)
									else: ax3.plot(float(gestart[exon])/resolution+(float(gestart[exon+1])/resolution-float(geend[exon])/resolution)/2+0.125,trackCount-gtrack+0.125, marker=arrow, color=icolor, markersize=1.25)					
												
					else:
							
						if len(genes[item])>5: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=genes[item][5])
						else: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=gcolor)
						ax3.add_patch(rect)			
								
				else: # simple plotting
										
					gstart = float(item.split('-')[0])/resolution
					gend = float(item.split('-')[1])/resolution
					gtrack = genes[item][0]
					
					if len(genes[item])>5: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=genes[item][5])
					else: rect = Rectangle((gstart,trackCount-gtrack), (gend-gstart), 0.25, color=gcolor)
					ax3.add_patch(rect)
									
 				
 				if plength <= 30 and geneLabels: # also consider the gene density
					
					#optimize the font size	
 				
 					gindex = nearest[gtrack].index(int(item.split('-')[1]))
 					upgene = nearest[gtrack][gindex-1]
 					if gindex < len(nearest[gtrack])-1: downgene = nearest[gtrack][gindex+1] 
 					else: downgene = upgene
 
					if plength <= 2 or plength < 1: plength=1
					elif plength <= 4: plength = 2
					else : plength/1.5
					gdist = min(abs(nearest[gtrack][gindex]-upgene),abs(nearest[gtrack][gindex]-downgene))
					if len(nearest[gtrack])==1: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=4.5/plength)
					elif float(gdist)/resolution >= 2 and len(genes[item][1])<=6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=4.5/plength)
					elif float(gdist)/resolution >= 2 and len(genes[item][1])>6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=3/plength)
					elif float(gdist)/resolution < 2 and float(gdist)/resolution > 1 and len(genes[item][1])<=6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=2.5/plength)
					elif float(gdist)/resolution < 2 and float(gdist)/resolution >1 and len(genes[item][1])>6: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=2/plength)
					elif float(gdist)/resolution <= 1 and float(gdist)/resolution >= 0.25: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=1.8)
					else: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=1)
				
					#if lblstndrd == 1:
 						#if len(genes[item][1])>5: ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=0.5)
 						#else : ax3.text(gstart, trackCount-gtrack+0.5, genes[item][1], fontsize=3)	
			
			ax3.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
			ax3.set_ylim(0,trackCount+1)
			plt.setp(ax3.get_yticklabels(), visible=False)
			if plotTadDomains and plotDomainTicks:
				ax3.set_xticks(tricks, minor=True)
				ax3.xaxis.grid(True,which='minor')
			
			if h_start > 0:
				for item in range(0,len(h_start)):
					ax3.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
			
			ax3.spines['right'].set_visible(False)
			ax3.spines['left'].set_visible(False)
			ax3.spines['top'].set_visible(False)
			ax3.tick_params(left="off")
			ax3.tick_params(right="off")
			ax3.xaxis.set_ticks_position('bottom')
						
			crowcounter+=2
			if numOfrows <= crowcounter and not randomBins: ax3.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))
			elif numOfrows <= crowcounter and randomBins: ax3.set_xlabel('Chromosome %s (Genomic Bins)' % (schr))
	
			
	# Pair-wise comparisons
	if compare and pair:
		
		fig.subplots_adjust(hspace=1.0)
		
		for peer1 in range(0,len(files)):
			
			prowcounter = rowcounter
			
			matrix1 = marray[peer1]
			
			if len(compareSm)>0 :
				slow = int(compareSm.split(',')[0])
				shigh = int(compareSm.split(',')[1])
				matrix1[(matrix1>=slow) & (matrix1<=shigh)]=1
			
			for peer2 in range(0,len(files)):	
				
				if peer1 != peer2:
					
					matrix2 = marray[peer2]
					
					if len(compareSm)>0 :
						slow = int(compareSm.split(',')[0])
						shigh = int(compareSm.split(',')[1])
						matrix2[(matrix2>=slow) & (matrix2<=shigh)]=1
					

					with np.errstate(divide='ignore',invalid='ignore'): matrix=matrix1/matrix2
					
					pax1 = plt.subplot2grid((numOfrows, 4*len(files)), (prowcounter, peer1*4), rowspan=4,colspan=4,sharex=ax1)
					with np.errstate(divide='ignore'): img=pax1.imshow(log2(matrix),cmap=plt.get_cmap("bwr"),origin="lower",interpolation="nearest",extent=(int(start or 1) - 0.5,\
																			  int(start or 1) + length - 0.5,int(start or 1) - 0.5,int(start or 1) + length - 0.5),aspect='auto')
		
					pax1.set_ylim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
					pax1.set_xlim(int(start or 1) - 0.5,int(start or 1) + length - 0.5)
					pax1.set_title(('log2(%s / %s)') % (names[peer1],names[peer2]))
					
					prowcounter+=4
					
					if len(peakFiles) > 0:
						origin_x,origin_y,radius,colors = read_peakFile(peakFiles[peer1],resolution,chromosome)
						for citem in range(0,len(origin_x)):
							if len(colors)==0: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor='black', linewidth=1, alpha=0.85)
							else: circle = Circle((origin_x[citem], origin_y[citem]), radius[citem], facecolor='none', edgecolor=colors[citem], linewidth=3, alpha=0.85)
							pax1.add_patch(circle)
				
					divider = make_axes_locatable(pax1)
					if len(compareEx)>0 : clow = int(compareEx.split(',')[0])*-1;chigh=int(compareEx.split(',')[1]);img.set_clim([clow,chigh])
					else: img.set_clim([-4,4])
				
					if wholeGenome : plt.setp(pax1.get_yticklabels(), visible=False)
					pax1.get_yaxis().set_label_coords(-0.125,0.5) 
					if plotTadDomains and plotDomainTicks:
						pax1.set_xticks(tricks, minor=True)
						pax1.xaxis.grid(True,which='minor',linewidth=2)
		
					if h_start > 0:
						for item in range(0,len(h_start)):
							pax1.axvspan(h_start[item], h_end[item], facecolor='g', alpha=0.10, linestyle='dashed')
		
					pax1.get_xaxis().set_label_coords(0.5,-0.125)
					
					if numOfrows <= prowcounter and not randomBins:
						pax1.set_xlabel('Chromosome %s Mb (resolution: %sKb)' % (schr , resolution/1000))	

					cax = divider.append_axes("bottom", size="2.5%", pad=0.1)
					cbar = plt.colorbar(img, cax=cax, ticks=MultipleLocator(2.0), format="%.1f",orientation='horizontal',extendfrac='auto',spacing='uniform')
					plt.setp(pax1.get_xticklabels(), visible=False)
		
					if not randomBins:
						ticks= pax1.get_xticks().tolist()
						if resolution*(end-start)<=500000: 
							for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,3)
						elif resolution*(end-start)<=1500000:
							for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,2)
						else: 
							for item in range(0,len(ticks)): ticks[item]=round(ticks[item]*resolution/1000000,1) 
						pax1.set_xticklabels(ticks)
						pax1.set_yticklabels(ticks)
						
					
	if len(oExtension) > 0 and oExtension in plt.gcf().canvas.get_supported_filetypes().keys(): extension='.'+oExtension
	elif 'JPEG' in plt.gcf().canvas.get_supported_filetypes_grouped().keys() or 'Joint Photographic Experts Group' in plt.gcf().canvas.get_supported_filetypes_grouped().keys(): extension='.jpeg'
	else : extension = '.png'
	
	warnings.simplefilter(action = "ignore", category = FutureWarning)
	
	print 'Plotting now!!'	
	if wholeGenome:
		if highResolution and dPixels==200:
			plt.savefig(output+'-WholeGenome-'+str(resolution/1000)+'K'+extension,dpi=200)
		elif dPixels !=200:
			plt.savefig(output+'-WholeGenome-'+str(resolution/1000)+'K'+extension,dpi=dPixels)
		else:
			plt.savefig(output+'-WholeGenome-'+str(resolution/1000)+'K'+extension)
	elif randomBins:
		if highResolution and dPixels==200:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').RandomBins'+extension,dpi=200)
		elif dPixels!=200:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').RandomBins'+extension,dpi=dPixels)
		else:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').RandomBins'+extension)
	else:
		if highResolution  and dPixels==200:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').'+str(resolution/1000)+'K'+extension,dpi=200)
		elif dPixels!=200:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').'+str(resolution/1000)+'K'+extension,dpi=dPixels)
		else:
			plt.savefig(output+'-'+chromosome+'.'+'ofBins('+str(start)+'-'+str(end)+').'+str(resolution/1000)+'K'+extension)

if __name__=='__main__':
	
	parser = argparse.ArgumentParser(usage='HiCPlotter.py -f file1 file2 ... -n name1 name2 ... -chr chr12 -o hES',add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
	
	group = parser.add_argument_group("Required Parameters")
	group.add_argument('-f','--files', nargs='+',help='',metavar='',required=True)
	group.add_argument('-n','--names', nargs='+',metavar='',required=True)
	group.add_argument('-chr', '--chromosome',default='',metavar='',required=True)
	group.add_argument('-o', '--output',default='',metavar='',required=True)
	
	group1 = parser.add_argument_group("Optional Parameters")
	group1.add_argument('-h', '--help', action="help")
	group1.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true",default=False)
	group1.add_argument('-d', '--dark',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-tri', '--tripleColumn',default=False,type=int,metavar='',help='default:0 - enable with 1')
	group1.add_argument('-bed', '--bedFile',default='',metavar='',help='')
	group1.add_argument('-hist', '--histograms', nargs='+',metavar='',default=[])
	group1.add_argument('-hl', '--histLabels', nargs='+',metavar='',default=[])
	group1.add_argument('-hm', '--histMax', nargs='+',metavar='',default=[])
	group1.add_argument('-fhist', '--fillHist', nargs='+',metavar='',default=[],help='(0:no, 1:yes)')
	group1.add_argument('-hc', '--histColors', nargs='+',metavar='',default=[])
	group1.add_argument('-si', '--superImpose',metavar='',type=int,default=False,help="default: 0 - enable with 1")
	group1.add_argument('-b', '--barPlots', nargs='+',metavar='',default=[])
	group1.add_argument('-bl', '--barLabels', nargs='+',metavar='',default=[])
	group1.add_argument('-bc', '--barColors', nargs='+',metavar='',default=[])
	group1.add_argument('-bm', '--barMax', nargs='+',metavar='',default=[])
	group1.add_argument('-t', '--tilePlots', nargs='+',metavar='',default=[])
	group1.add_argument('-tl', '--tileLabels', nargs='+',metavar='',default=[])
	group1.add_argument('-tc', '--tileColors', nargs='+',metavar='',default=[])
	group1.add_argument('-tt', '--tileText',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-a', '--arcPlots', nargs='+',metavar='',default=[])
	group1.add_argument('-al', '--arcLabels', nargs='+',metavar='',default=[])
	group1.add_argument('-ac', '--arcColors', nargs='+',metavar='',default=[])
	group1.add_argument('-high', '--highlights',default=0,type=int,metavar='',help='default:0 - enable with 1')
	group1.add_argument('-hf', '--highFile',default='',metavar='',help='')
	group1.add_argument('-peak', '--peakFiles', nargs='+',metavar='',default=[])
	group1.add_argument('-g', '--plotGenes', nargs='+',metavar='',default='')
	group1.add_argument('-gl', '--geneLabels',type=int,default=True,metavar='',help="default: 1 - disable with 0")
	group1.add_argument('-ep', '--epiLogos',metavar='',default='')
	group1.add_argument('-ext', '--oExtension',default='',metavar='')
	group1.add_argument('-spi', '--spine',metavar='',type=int,default=False,help="default: 0 - enable with 1")
	group1.add_argument('-im', '--imputed',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-s', '--start',type=int,default=0,metavar='',help="default: 0")
	group1.add_argument('-e', '--end',type=int,default=0,metavar='',help="default: matrix end")
	group1.add_argument('-r', '--resolution',type=int,default=100000,metavar='',help="default: 100000")
	group1.add_argument('-rb', '--randomBins',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-wg', '--wholeGenome',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-w', '--window',type=int,default=5,metavar='',help="default: 5")
	group1.add_argument('-fh', '--fileHeader',type=int,default=0,metavar='',help="default: 0")
	group1.add_argument('-ff', '--fileFooter',type=int,default=0,metavar='',help="default: 0")
	group1.add_argument('-tr', '--tadRange',type=int,default=8,metavar='',help="default: 8")
	group1.add_argument('-hmc', '--heatmapColor',type=int,default=3,metavar='',help="Colors for heatmap: Greys(0), Reds(1), YellowToBlue(2), YellowToRed(3-default), Hot(4), BlueToRed(5)")
	group1.add_argument('-sn', '--smoothNoise',type=float,default=0.5,metavar='',help="default: 0.5")
	group1.add_argument('-mm', '--matrixMax',type=int,default=10,metavar='',help="default: 0")
	group1.add_argument('-dpi', '--dPixels',type=int,default=200,metavar='',help="default: 0")
	group1.add_argument('-c', '--compare',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-ce', '--compareEx',default='',metavar='',help='')
	group1.add_argument('-cs', '--compareSm',default='',metavar='',help='')
	group1.add_argument('-p', '--pair',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-cn', '--cleanNANs',type=int,default=True,metavar='',help="default: 1 - disable with 0")
	group1.add_argument('-hR', '--highResolution',type=int,default=True,metavar='',help="default: 1 - disable with 0")
	group1.add_argument('-pi', '--plotInsulation',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-dc', '--domColors', nargs='+',metavar='',default=[])
	group1.add_argument('-ptr', '--plotTriangular',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-trh', '--triangularHeight',type=int,default=False,metavar='',help="default: (end-start)/5")
	group1.add_argument('-ptd', '--plotTadDomains',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-pdt', '--plotDomainTicks',type=int,default=True,metavar='',help="default: 1 - enable with 0")
	group1.add_argument('-pcd', '--plotCustomDomains',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-pdb', '--plotDomainsAsBars',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-pcdf', '--customDomainsFile',nargs='+',metavar='',default=[])
	group1.add_argument('-pptd', '--plotPublishedTadDomains',type=int,default=False,metavar='',help="default: 0 - enable with 1")
	group1.add_argument('-ptdo', '--publishedTadDomainOrganism',type=int,default=True,metavar='',help="human(default): 1 - mouse: 0")
	
	args = vars(parser.parse_args())
	
	if len(args['files']) != len(args['names']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and names'
		raise SystemExit
	if len(args['histograms'])>0 and len(args['histograms'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and BedGraphs'
		raise SystemExit
	if len(args['histLabels'])>0 and len(args['histLabels'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and BedGraph Labels'
		raise SystemExit
	if len(args['histograms'])>0:
		if len(args['histograms'][0].split(','))>2 and args['superImpose']:
			print >>sys.stderr, 'Upps!! Please use super impose (-si) only with 2 histogram files per condition'
			raise SystemExit
	if len(args['histograms'])+len(args['histLabels'])+len(args['fillHist'])>0 and len(args['histLabels'])!=len(args['histograms']) and len(args['histLabels'])!=len(args['fillHist']):
		print >>sys.stderr, 'Upps!! Please provide equal number of BedGraphs, BedGraph Labels and FillUnders (0:no, 1:yes)'
		raise SystemExit
	if len(args['barPlots'])>0 and len(args['barPlots'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and bar plots'
		raise SystemExit
	if len(args['barLabels'])>0 and len(args['barLabels'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and bar plot Labels'
		raise SystemExit
	if len(args['barPlots'])+len(args['barLabels'])>0 and len(args['barPlots'])!=len(args['barLabels']):
		print >>sys.stderr, 'Upps!! Please provide equal number of bar plot and bar plot Labels'
		raise SystemExit
	if len(args['barColors'])>0 and len(args['barPlots'])!=len(args['barColors']):
		print >>sys.stderr, 'Upps!! Please provide equal number of bar plot and bar plot colors'
		raise SystemExit
	if len(args['tilePlots'])>0 and len(args['tilePlots'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and tile plots'
		raise SystemExit
	if len(args['tileLabels'])>0 and len(args['tileLabels'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and tile plot Labels'
		raise SystemExit
	if len(args['tilePlots'])+len(args['tileLabels'])>0 and len(args['tilePlots'])!=len(args['tileLabels']):
		print >>sys.stderr, 'Upps!! Please provide equal number of tile plot and tile plot Labels'
		raise SystemExit
	if len(args['tileColors'])>0 and len(args['tilePlots'])!=len(args['tileColors']):
		print >>sys.stderr, 'Upps!! Please provide equal number of tile plot and tile plot colors'
		raise SystemExit
	if len(args['arcPlots'])>0 and len(args['arcPlots'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and arc plots'
		raise SystemExit
	if len(args['arcLabels'])>0 and len(args['arcLabels'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and arc plot Labels'
		raise SystemExit
	if len(args['arcPlots'])+len(args['arcLabels'])>0 and len(args['arcPlots'])!=len(args['arcLabels']):
		print >>sys.stderr, 'Upps!! Please provide equal number of arc plot and arc plot Labels'
		raise SystemExit
	if len(args['arcColors'])>0 and len(args['arcPlots'])!=len(args['arcColors']):
		print >>sys.stderr, 'Upps!! Please provide equal number of arc plot and arc plot colors'
		raise SystemExit
	if args['plotCustomDomains'] and len(args['customDomainsFile'])==0:
		print >>sys.stderr, 'Upps!! Please provide a bedGraph file for custom domains'
		raise SystemExit
	if args['plotCustomDomains'] and len(args['customDomainsFile'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and custom domains'
		raise SystemExit
	if args['start'] < 0 or args['end'] < 0 or args['end'] - args['start'] < 0:
		print >>sys.stderr, 'Upps!! Start and end should be positive and end bigger than start'
		raise SystemExit
	if len(args['peakFiles'])>0 and len(args['peakFiles'])!=len(args['files']):
		print >>sys.stderr, 'Upps!! Please provide equal number of HiC matrix and peak files'
		raise SystemExit
		
	if args['verbose']:
		logging.basicConfig(filename=args['output']+'.log',level=logging.DEBUG,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
		logging.info('You are using HiCPlotter version:%s',version)
		logging.info('Using arguments: %s',args)
		logging.info('\n#################################\n')
	
	HiCplotter(**args)