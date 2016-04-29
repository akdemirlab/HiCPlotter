import sys

def bedToSparse():
	
	'''
    generates triple sparse output from bed-formatted matrix files

    parameters:
	
	genome: a tab delimited file for chromosome sizes. [CHR Length] chr1	249250621
	windowsize: bin size of the bed file.
	bed : bed-formatted matrix file. Example line: chr10	355155_356426	483838_485372	1
	
	creates:
	bed file for the bins
	matrix file formatted as triple-sparse matrix (HiC-Pro output)
    
    Special thanks to Charles Dietz for testing this script!
    
    '''

	
	genome=sys.argv[1]
	windowsize=int(sys.argv[2])
	bed=sys.argv[3]
	
	writer = open(bed+'_'+str(windowsize/1000)+'Kb_ord.bed','w')
	counter = 1
	locs = {}
	
	fone=open(genome,'r')
	for line in fone.xreadlines():
		l=line.strip().split()
		locs[l[0]]=[]
		start=0
		while start < int(l[1]):
			locs[l[0]].append(counter)
			writer.write(l[0]+'\t'+str(start)+'\t'+str(start+windowsize)+'\t'+str(counter)+'\n')
			start = start+windowsize
			counter +=1
	
	
	writer = open(bed+'_'+str(windowsize/1000)+'Kb.matrix','w')
	fone=open(bed,'r')
	for line in fone.xreadlines():
		tabs=line.strip().split()
		anchor1 = int(tabs[1].split("_")[0])/windowsize
		anchor2 = int(tabs[2].split("_")[1])/windowsize
		writer.write(str(locs[tabs[0]][anchor1])+'\t'+str(locs[tabs[0]][anchor2])+'\t'+tabs[3]+'\n')
		
if __name__ == "__main__":
    bedToSparse()
