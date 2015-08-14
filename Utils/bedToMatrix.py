from array import *
import math
import sys

def bedToMatrix(): 
	fbed=sys.argv[1]
	chromosome = {}
	max_freq = 0
	last_region = 0
	
	fone=open(fbed,'r')
	for line in fone.readlines():
		tabs=line.strip().split("\t")
		resolution = int(tabs[1].split("_")[1])-int(tabs[1].split("_")[0])
		start = int(tabs[1].split("_")[0])/resolution
		interact = int(tabs[2].split("_")[0])/resolution
		if last_region < int(tabs[2].split("_")[1])/resolution: last_region = int(tabs[2].split("_")[1])/resolution
		freq = float(tabs[3])
		if freq > max_freq: max_freq = freq
		if start not in chromosome.keys():
			chromosome[start]={}
			chromosome[start][interact]=freq
		else:
			chromosome[start][interact]=freq			
	fone.close()
	
	writer = open(fbed+'_matrix'+'.txt','w')
	for i in range(0,last_region):
		if i in chromosome.keys():
			for x in range(0,last_region):
				if x == i:
					if x not in chromosome[i].keys():
						writer.write(str(max_freq)+'\t')
					elif x in chromosome.keys() and i in chromosome[x].keys():
						writer.write(str(chromosome[x][i])+'\t')
					elif x in chromosome[i].keys():
						writer.write(str(chromosome[i][x])+'\t')
				else:
					if x in chromosome[i].keys():
						writer.write(str(chromosome[i][x])+'\t')
					elif x in chromosome.keys() and i in chromosome[x].keys():
						writer.write(str(chromosome[x][i])+'\t')
					else:
						writer.write('0\t')
		else:
			for x in range(0,last_region):
				if x == i:
					writer.write(str(max_freq)+'\t')
				else:
					writer.write('0\t')
		writer.write('\n')
	
				
if __name__ == "__main__":
    bedToMatrix()