#HiCPlotter: Integrating HiC data with genomic datasets

# About HiCPlotter:

HiCPlotter is a Python data visualization tool for integrating different data types with interaction matrixes. For more on 5C or HiC, please check:  <a class="reference external" href="http://www.nature.com/nrg/journal/v14/n6/full/nrg3454.html"> Dekker et al. 2013</a>.

_HiCPlotter is designed by Kadir Caner Akdemir (kcakedemir at mdanderson dot org) in Lynda Chin's Lab at the University of Texas MD Anderson Cancer Center, Houston, TX, USA._

# Requirements:

Python 2.7.*

[Numpy, Scipy, Matplotlib](http://www.scipy.org/).

* Please note: scipy, numpy and matplotlib modules should be installed and updated to current version. If you get error(s) related to one or more of these modules, check this [solution](https://github.com/kcakdemir/HiCPlotter/issues/1).

_HiCPlotter is tested on Mac OS (Mountain Lion and Yosemite) and Linux (RedHat 4.1.2-44 and 5.5-Final) systems._

_HiCPlotter is purposefully designed with the least amount of dependencies to make it easily applicable._


# Arguments:

	Required parameters:

    files: a list of filenames to be plotted.
    name: a list of labels for the experiment.
    chr: chromosome to be plotted.
    output: prefix for the output file.
    
    Optional parameters:
    
    verbose: print version and arguments into a file
    histograms: a list of filenames to be plotted as histogram.
    histLabels: a list of labels for the histograms.
    fillHist: a list whether each histogram will be filled (1) or not (0:default).
    histMax : a list of integer for maximum values of histograms.
    start: retain after x-th bin (0:default).
    end: continues until x-th bin (default: length of the matrix).
    resolution: resolution of the bins (default: 100000).
    tilePlots: a list of filenames to be plotted as tile plots.
    tileLabels: a list of labels for the tile plots.
    tileColors: a list of hexadecimal numbers for coloring the tile plots.
    tileText: an integer whether text will be displayed above tiles (0:default) or not (1).
    arcPlots: a list of filenames to be plotted as arc plots.
    arcLabels: a list of labels for the arc plots.
    arcColors: a list of hexadecimal numbers for coloring the arc plots.
    highlights: an integer for enabling highlights on the plot (0:default), enable(1). 
    highFile: a file name for a bed file to highlight selected intervals.
    peakFiles : a list of filenames to be plotted on the matrix.
    window: an integer of distance to calculate insulation score.
    tadRange: an integer of window to calculate local minima for TAD calls.
    fileHeader: an integer for how many lines should be ignored in the matrix file (1:default).
    fileFooter: an integer for how many lines should be skipped at the end of the matrix file (0:default).
    smoothNoise: a floating-point number to clean noise in the data.
    heatmapColor: an integer for choosing heatmap color codes: Greys(0), Reds(1), YellowToBlue(2), YellowToRed(3-default), Hot(4), BlueToRed(5).
    cleanNANs: an integer for replacing NaNs in the matrix with zeros (1:default) or not (0).
    plotTriangular: an integer for plotting rotated half matrix (1:default) or not (0).
    plotTadDomains: an integer for plotting TADs identified by HiCPlotter (1) or not (0:default).
    plotPublishedTadDomins: an integer for plotting TADs from Dixon et, al. 2012 (1:default) or not (0).
    highResolution: an integer whether plotting high resolution (1:default) or not (0).
    plotInsulation: an integer for plotting insulation scores (0:default) or plot (1).
    randomBins: an integer for plotting random resolution data (1:default) or not (0).
    wholeGenome: an integer for plotting whole genome interactions (1:default) or not (0).
    plotCustomDomains: a list of file names to be plotted beneath the matrix.
    publishedTadDomainOrganism : an integer for plotting human (1:default) or mouse (0) TADs from Dixon et, al. 2012.
    customDomainsFile: a list of filenames to be plotted as TADs for each experiments.

# Input Files

## Hi-C/5C matrix data

For visualizing Hi-C or 5C data, HiCPlotter requires a matrix file (by default first line is ignored).

Matrix files are plotted as their log2 values and color legend is put below the plot.

	Bin1	Bin2	Bin3	Bin4	Bin5	Bin6
	7.85957	4.80329	11.4766	9.57416	4.5288	8.55022
	8.61621	4.98956	2.35654	5.69483	11.1187	10.1322
	4.06803	4.07801	7.98047	2.59144	6.3851	7.74306
	4.52869	2.70624	8.94544	4.29185	8.29491	8.38257
	2.91472	3.84658	1.56752	4.48515	7.4955	8.77461
	3.08096	2.96487	7.23623	2.33142	3.08529	5.5379
	3.12141	3.06905	4.97247	2.39298	5.03621	7.22344
	3.4037	2.26455	1.48176	1.41958	3.40252	7.7027
	3.8696	1.41425	7.68872	2.21027	5.06846	3.20063

If your file (example is modified from [GSM873926](www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM873926)) has several header lines, use -fh parameter (this particular case -fh 5).
	
	# resultName	mESCs-female-PGK-day2-Replicate1			
	#				
	#				
	#
	REV_2|mm9|chrX:98831149-98834145	REV_4|mm9|chrX:98837507-98840771	REV_6|mm9|chrX:98841228-98843248	REV_12|mm9|chrX:98855723-98862021
	936.010581657246					743.499513378904					241.956223702097					23.2451328286973
	69.8831429744098					513.412096905019					747.143877424081					7.1902317648089

## BedGraph

For visualizing any type of genomic data, HiCPlotter uses bedGraph format.
	
	
	chromA  chromStartA  chromEndA  dataValueA color     text
	chr1	10000		 10500		10.0	   250,13,27 Polycomb

4th column should be a floating number for histograms.

5th column should be an rgb color for tile or arc plots. 

6th column should be a string for tile plots. _Columns after 6th will be ignored._

## Peak File

For annotating the interaction matrix, HiCPlotter requires the following format.
	
	chromA	chromStartA  chromEndA	chromB	chromStartB chromEndB	color (optional)
	chr10	100180000	100190000	chr10	100410000	100420000	0,255,255
	chr10	101600000	101610000	chr10	101800000	101810000	0,255,255
	chr10	102100000	102105000	chr10	102190000	102195000	0,255,255

# Usage:

	python HiCPlotter.py -f file1 file2 -n name1 name2 -chr chrX -o output

_**You can generate all of the following examples by using "sh testRun.sh", all of required files are in the data folder.**_

## Plotting whole chromosome
_Please note: -fh is set to 0 as the input matrix doesn't have a header line._

_Hi-C and TADs data taken from:_ [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517)

	python HiCPlotter.py -f data/HiC/Human/hES-nij.chr21.2 -n hES -chr chr21 -r 40000 -o default1 -fh 0

<figure>
  <figcaption align="middle">**Chromosome interactions at chromosome 21 in human ES cells**</figcaption>
  <img src="examplePlots/default1-chr21.ofBins(0-1174).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Focusing to a region within a chromosome
_Start and end locations can be specified as bin numbers with -s and -e parameters._

_Color of triangles specify interaction frequency in a given TAD._

_TADs identified by Dixon et al. can be plotted with -pptd parameter._	
	
	python HiCPlotter.py -f data/HiC/Human/hES-nij.chr21.2 -n hES -chr chr21 -r 40000 -o default2 -ptd 1 -pptd 1 -s 600 -e 900 -fh 0 -w 8 -tr 10 -pi 1

<figure>
  <figcaption align="middle">**Chromosome interactions at chromosome 21 in human ES cells**</figcaption>
  <img src="examplePlots/default2-chr21.ofBins(600-900).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Visualization of multiple Hi-C datasets

_Hi-C data taken from:_ [Zuin et, al. PNAS 2014](http://www.pnas.org/content/111/3/996.long)

_Color code of the heatmaps can be changed with -hmc parameter_
	
 	python HiCPlotter.py -f data/HiC/Human/GSM1081526_TEV_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081528_HRV_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081530_CTRL_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081533_CTCF_r2_cis.index.chr6.txt_matrix.txt -n WT RAD21-Depleted siControl CTCF-Depleted -chr chr6 -r 40000 -fh 0 -pi 0 -sn 0.35 -o Rad21.CTCF -s 2800 -e 2950 -hmc 5

<figure>
  <figcaption align="middle">**Chromosome interactions for wild type, RAD21-depleted, CTCF-depleted cells.**</figcaption>
  <img src="examplePlots/Rad21.CTCF-chr6.ofBins(2800-2950).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>


# Example cases with publicly available datasets:

## Visualization of ChIP-Seq and 4C data as histograms

_Multiple histograms for the same matrix should be seperated by comma (true for hist labels and fill histogram parameters)._

_Data taken from:_ 4C : [Noordermer et, al. Elife 2014](http://elifesciences.org/content/3/e02557), Hi-C and TADs : [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517) and CTCF : [Stadler et, al. Nature 2011](http://www.nature.com/nature/journal/v480/n7378/full/nature10716.html)

	python HiCPlotter.py -f data/HiC/Mouse/mES.chr2 -n mES -chr chr2 -r 40000 -o HoxD -hist data/HiC/Mouse/GSM1334415_4C_Mouse_EScells_Hoxd4_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334440_4C_Mouse_E9.5TB_Hoxd4_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334412_4C_Mouse_EScells_Hoxd13_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334437_4C_Mouse_E9.5TB_Hoxd13_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM747534_ChIPseq_CTCF_ES_rep1.chr2.bedGraph -hl Hoxd4-ES,Hoxd4-Tail,Hoxd13-ES,Hoxd13-Tail,CTCF-ES -s 1830 -e 1880 -fh 0 -pi 0 -pcd 1 -pcdf data/mES_domains_mm9.bed -fhist 1,1,1,1,0 -hm 2000,2000,2000,2000,50

<figure>
  <figcaption align="middle">**3D compartments in the HoxD cluster in mouse ES cells**</figcaption>
  <img src="examplePlots/HoxD-chr2.ofBins(1830-1880).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Visualization of ChIP-Seq and RAP-Seq data as histograms

_Data taken from:_ RAP-seq : [Engreitz et al. Science 2014](http://www.sciencemag.org/content/341/6147/1237973.long), Hi-C : [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517) and H3K27me3 : [Mouse ENCODE Project](http://www.mouseencode.org/)

_Rotated matrix can be removed with -ptr 0 parameter_

 	python HiCPlotter.py -f data/HiC/Mouse/mES.chrX -n mES -r 40000 -chr chrX -o RAP -fh 0 -hist data/HiC/Mouse/GSE46918_pSM33-0hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-1hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-2hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-3hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-6hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/wgEncodeLicrHistoneEsb4H3k27me3ME0C57bl6StdSig.chrX.bedGraph -hl Xist_0h,Xist_1h,Xist_2h,Xist_3h,Xist_6h,H3K27me3_0h -pi 0 -ptr 0 -fhist 0,1,1,1,1,0 -hmc 4 -sn 0

<figure>
  <figcaption align="middle">**Xist spreading during initiation of X-chromosome inactivation**</figcaption>
  <img src="examplePlots/RAP-chrX.ofBins(0-4167).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Visualization of ChIP-Seq as histograms, ChIA-Pet as arcs and Polycomb domains as tiles

_Arc plots require a bedGraph file (-a file1), color can be specied as a hexadecimal number (-ac B4B4B4) or for each arc by specified RGB colors in bedGraph file._

_Data taken from:_ SMC ChIA-Pet and Polycomb Domains: [Dowen et, al. Cell 2014](http://www.sciencedirect.com/science/article/pii/S0092867414011799), Hi-C and TADs : [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517) and H3K27me3 : [Mouse ENCODE Project](http://www.mouseencode.org/)
	
	python s.py -f data/HiC/Mouse/mES.chr3 -n mES -chr chr3 -o Bhlhe22 -r 40000 -s 400 -e 500 -a data/HiC/Mouse/mESC_SMC_ChIPPet.bed -al SMC -hist data/HiC/Mouse/GSM747534_chr3.bedGraph,data/HiC/Mouse/wgEncodeLicrHistoneEsb4H3k27me3ME0C57bl6StdSig.chr3.bedGraph -hl CTCF,H3K27me3 -pi 0 -ptr 0 -t data/HiC/Mouse/mm9_Polycomb_domains.bed -tl Polycomb -tc 00CCFF -ac B4B4B4 -fh 0

<figure>
  <figcaption align="middle">**Bhlhe22 locus in mouse ES cells**</figcaption>
  <img src="examplePlots/Bhlhe22-chr3.ofBins(400-475).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Visualization of 4C data as histograms and Enhancers as tiles with text

_If bedGraph file for tile plotting contains text in 6th column, features can be plotted above tiles with -tt parameter._

_Data taken from:_ 4C : [Lonfat et, al. Science 2014](http://www.sciencemag.org/content/346/6212/1004.long) (Primer locations in Lonfat et, al. are extended to make them more visible in the plot), Hi-C and TADs : [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517)


	python HiCPlotter.py -f data/HiC/Mouse/mES.chr6 -n mES -chr chr6 -r 40000 -o Digit.vs.GT -s 1295 -e 1338 -hist data/HiC/Mouse/GSM1524258_segToFrag_4C_Digits_WT_E12-5_HoxA13_smoothed_11FragsPerWin.bedGraph,data/HiC/Mouse/GSM1524259_segToFrag_4C_GT_WT_E15-5_HoxA13_smoothed_11FragsPerWin.bedGraph -hl Digits,GT -fhist 1,1 -fh 0 -pi 0 -hm 1500,1500 -pcd 1 -pcdf data/mES_domains_mm9.bed -sn 0.4 -t data/HiC/Mouse/primers.bedGraph -tl Enhancers -tt 1

<figure>
  <figcaption align="middle">**Interaction profiles of Hoxa13 gene with tissue specific enhancers**</figcaption>
  <img src="examplePlots/Digit.vs.GT-chr6.ofBins(1295-1338).40K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Highlighting selected loci on the plot

_Highlights on the plots can be drawn with -high 1 and passing a bed file name to -hf parameter._

_Data taken from:_ Hi-C and Arrow domains  : [Rao et, al. Cell 2014](http://www.cell.com/cell/abstract/S0092-8674\(14\)01497-4), ChIA-Pet [Heidari et, al. Genome Research 2014](http://genome.cshlp.org/content/early/2014/09/15/gr.176586.114), TADs : [Dixon et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11082.html?WT.ec_id=NATURE-20120517), ChIP-Seq and Repli-Seq data : [Encode Project](http://www.genome.gov/encode/)


	python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -chr chr10 -r 25000 -pi 0 -fh 0 -o ChIA -a data/HiC/Human/GM12878.Rad21.bed data/HiC/Human/K562.Rad21.bed -al ChIA-PET ChIA-PET -s 3000 -e 3500 -pcd 1 -pcdf data/HiC/Human/GM12878_Arrowhead_domainlist.bed data/HiC/Human/K562_Arrowhead_domainlist.bed -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -hl DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq -fhist 0,0,1 0,0,1 -pptd 1 -high 1 -hf data/HiC/Human/highlight.bed

<figure>
  <figcaption align="middle">**Interaction profiles in chromosome 10 of GM12878 and K562 cells.**</figcaption>
  <img src="examplePlots/High-chr10.ofBins(3000-3500).25K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Annotating the interaction matrix

_Annotations on the matrix can be drawn with -peak parameter.Input file should contain at least six columns._

_Data taken from:_ Hi-C and HiCCUP peaks  : [Rao et, al. Cell 2014](http://www.cell.com/cell/abstract/S0092-8674\(14\)01497-4)


 	python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/KBM7-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HMEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 KBM7 K562 HUVEC IMR90 HMEC NHEK -r 25000 -pi 0 -fh 0 -o Loops -chr chr10 -peak data/HiC/Human/GSE63525_GM12878_replicate_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_KBM7_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_K562_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_HUVEC_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_IMR90_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_HMEC_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_NHEK_HiCCUPS_looplist.bed -s 3600 -e 3675  -ptr 0

<figure>
  <figcaption align="middle">**Significantly interaction loci (peaks) in several cell lines.**</figcaption>
  <img src="examplePlots/Loops-chr10.ofBins(3600-3675).25K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Whole Genome Plotting

_Whole genome plotting can be activated by -wg parameter (Please note: currently only matrixes can be plotted with this option)._

_Data taken from:_ Hi-C : [Seitan et, al. Genome Research 2014](http://genome.cshlp.org/cgi/pmidlookup?view=long&pmid=24002784)
	
 	python HiCPlotter.py -f data/HiC/Human/GSM1184323-HiCMYZ-Tcell-Rad21WT-R1.mm9.NA.L-1400000-wDiag-noSS-iced.2.matrix data/HiC/Human/GSM1184321-HiCMYZ-Tcell-Rad21KO-R1.mm9.NA.L-1400000-wDiag-noSS-iced.2.matrix -n Tcell_WT Tcell_Rad21KO -chr Genome -r 1400000 -o Tcell -pi 0 -ptr 0 -wg 1 -hmc 5 -fh 4

<figure>
  <figcaption align="middle">**T-cell whole genome interaction data in wild type and Rad21 knock-out cells**</figcaption>
  <img src="examplePlots/Tcell-WholeGenome-1400K.jpeg" alt="Example plot from HiCPlotter">
</figure>

## 5C data visualization

_Random binned 5C data plotting can be activated by -rb parameter (Please note: currently only matrixes and triangular plots can be plotted with this option)._

_Data taken from:_ 5C data [Nora et, al. Nature 2012](http://www.nature.com/nature/journal/v485/n7398/full/nature11049.html)

 	python HiCPlotter.py -f data/5C/GSM873926_mESCs-female-PGK12.1-day2-Replicate1.txt data/5C/GSM873932_femaleXO-mESCs-DXTX-replicate-1.matrix.txt data/5C/GSM873924_female-MEFs-replicate-1.matrix.txt -n mESC mESC_XO MEF -fh 8 -chr chrX -o 5C -sn 2 -pi 0 -rb 1 -e 300 -hmc 5

<figure>
  <figcaption align="middle">**Xist locus in mouse ES, mouse ES Xist deletion and MEF cells**</figcaption>
  <img src="examplePlots/5C-chrX.ofBins(0-300).RandomBins.jpeg" alt="Example plot from HiCPlotter">
</figure>

# Figures

## Basic usage

	python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 HUVEC NHEK IMR90 -chr chr10 -r 25000 -s 3000 -e 3250 -o Figure1 -fh 0

<figure>
  <figcaption align="middle">**Basic Usage of HiCPlotter**</figcaption>
  <img src="examplePlots/Figure1.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Adding tracks

	python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 HUVEC NHEK IMR90 -chr chr10 -r 25000 -s 3000 -e 3500 -o Figure2 -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseHuvecRawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneHuvecCtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseNhekRawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneNhekCtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqNhekWaveSignalRep1.bedGraph data/HiC/Human/wgEncodeOpenChromDnaseImr90BaseOverlapSignal.chr10.bedGraph,data/HiC/Human/wgEncodeSydhTfbsImr90CtcfbIggrabSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqImr90WaveSignalRep1.bedGraph -fh 0 -fhist 0,0,1 0,0,1 0,0,1 0,0,1 0,0,1 -hl DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq -hm 400,600,100 400,600,100 400,600,100 400,600,100 400,600,100 -pcd 1 -pcdf data/HiC/Human/GM12878_Arrowhead_domainlist.bed data/HiC/Human/K562_Arrowhead_domainlist.bed data/HiC/Human/HUVEC_Arrowhead_domainlist.bed data/HiC/Human/NHEK_Arrowhead_domainlist.bed data/HiC/Human/IMR90_Arrowhead_domainlist.txt -t data/HiC/Human/GM12878_18_core_K27ac_dense2.bed data/HiC/Human/K562_18_core_K27ac_dense2.bed data/HiC/Human/HUVEC_18_core_K27ac_dense2.bed data/HiC/Human/NHEK_18_core_K27ac_dense2.bed data/HiC/Human/IMR90_18_core_K27ac_dense2.bed -tl ChromHMM ChromHMM ChromHMM ChromHMM ChromHMM -pptd 1 -high 1 -hf data/HiC/Human/fig2.bed 

<figure>
  <figcaption align="middle">**Adding tracks**</figcaption>
  <img src="examplePlots/Figure2.jpeg" alt="Example plot from HiCPlotter">
</figure>

## Cohesin ChIA-PET interactions coincide with early replication sites

	python HiCPlotter.py -f data/HiC/Human/GM12878-chr15_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr15_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -chr chr15 -r 25000 -s 1800 -e 2250 -o Figure3 -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep1.chr15.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr15.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep1.chr15.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr15.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -fh 0 -fhist 0,0,1 0,0,1 -hl DNase,CTCF,RepliSeq DNase,CTCF,RepliSeq -hm 400,400,100 400,400,100 -t data/HiC/Human/GM12878_Enhancer.bed,data/HiC/Human/GM12878_Txn.bed,data/HiC/Human/GM12878_Het.bed data/HiC/Human/K562_Enhancer.bed,data/HiC/Human/K562_Txn.bed,data/HiC/Human/K562_Het.bed -tl Enhancer,Transcribed,Heterochromatin Enhancer,Transcribed,Heterochromatin -a data/HiC/Human/GM12878.Rad21.bed data/HiC/Human/K562.Rad21.bed -al RAD21 RAD21 -ptr 0 -high 1 -hf data/HiC/Human/fig3.bed

<figure>
  <figcaption align="middle">**Cohesin ChIA-PET interactions coincide with early replication sites**</figcaption>
  <img src="examplePlots/Figure3.jpeg" alt="Example plot from HiCPlotter">
</figure>

## A lincRNA locus exhibits active chromatin formation in K562 cells

	python HiCPlotter.py -f data/HiC/Human/GM-chr19_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr19_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -r 25000 -chr chr19 -hist data/HiC/Human/GM12878.DNAse.chr19.2.bedGraph,data/HiC/Human/GM12878.RnaSeq.chr19.2.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/K562.DNAse.chr19.2.bedGraph,data/HiC/Human/K562.RnaSeq.chr19.2.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -hl DNAse,RNASeq,RepliSeq DNAse,RNASeq,RepliSeq -t data/HiC/Human/GM12878_TSS+Trx.2.bed data/HiC/Human/K562_TSS+Trx.2.bed -tl ChromHMM ChromHMM -high 1 -hf data/HiC/Human/region.bed -o Figure4 -s 1100 -e 1302 -hm 300,300,100 300,300,100 -fh 0 -fhist 0,0,1 0,0,1 -ptr 0

<figure>
  <figcaption align="middle">**A lincRNA locus exhibits active chromatin formation in K562 cells**</figcaption>
  <img src="examplePlots/Figure4.jpeg" alt="Example plot from HiCPlotter">
</figure>

# Tips

If your data contains several columns before data matrix, from command line you could use: 
	
	cut -f N- matrix > new_matrix (where N is the ith column data values start)

If any of the imported packages are missing in your python system, try to comment out those lines. For example:

	Original  :     from scipy.signal import argrelextrema (line 20)
	Try this  :     #from scipy.signal import argrelextrema (line 20). Use HiCPlotter with the -pi 0 and -ptd 0

If you like to run HiCPlotter in verbose mode, please use -v parameter which will create a log file with which parameters the program ran.

# Help:

If you encounter any problems, please contact with - Kadir Akdemir (kcakedemir at mdanderson dot org) or Dr. Lynda Chin (lchin at mdanderson dot org).

# Acknowledgements

Thanks to Lynda Chin for her leadership, management and support.

Thanks to Zeynep Coban-Akdemir, Ian Watson, Kunal Rai and Samir Amin for their insightful comments.

We are grateful to many researchers cited above for providing their data in publicly available and easy-to-use format.
