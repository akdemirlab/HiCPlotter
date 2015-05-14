#!/bin/bash

find ./data/ -name "*.gz" -exec gunzip -v {} \;

echo -e "\n#### Basic examples #### \n"

echo -e "\nPlotting whole chromosome"
python HiCPlotter.py -f data/HiC/Human/hES-nij.chr21.2 -n hES -chr chr21 -r 40000 -o default1 -fh 0
 
echo -e "\nFocusing to a region within a chromosome"
python HiCPlotter.py -f data/HiC/Human/hES-nij.chr21.2 -n hES -chr chr21 -r 40000 -o default2 -ptd 1 -pptd 1 -s 600 -e 900 -fh 0 -w 8 -tr 10 -pi 1

echo -e "\nVisualization of multiple Hi-C datasets"
python HiCPlotter.py -f data/HiC/Human/GSM1081526_TEV_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081528_HRV_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081530_CTRL_r1_cis.index.chr6.txt_matrix.txt data/HiC/Human/GSM1081533_CTCF_r2_cis.index.chr6.txt_matrix.txt -n WT RAD21-Depleted siControl CTCF-Depleted -chr chr6 -r 40000 -fh 0 -pi 0 -sn 0.35 -o Rad21.CTCF -s 2800 -e 2950 -hmc 5

echo -e "\n#### Examples from publicly available datasets ####"

echo -e "\nVisualization of ChIP-Seq and 4C data as histograms"
python HiCPlotter.py -f data/HiC/Mouse/mES.chr2 -n mES -chr chr2 -r 40000 -o HoxD -hist data/HiC/Mouse/GSM1334415_4C_Mouse_EScells_Hoxd4_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334440_4C_Mouse_E9.5TB_Hoxd4_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334412_4C_Mouse_EScells_Hoxd13_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM1334437_4C_Mouse_E9.5TB_Hoxd13_smoothed_11windows.bedGraph,data/HiC/Mouse/GSM747534_ChIPseq_CTCF_ES_rep1.chr2.bedGraph -hl Hoxd4-ES,Hoxd4-Tail,Hoxd13-ES,Hoxd13-Tail,CTCF-ES -s 1830 -e 1880 -fh 0 -pi 0 -pcd 1 -pcdf data/mES_domains_mm9.bed -fhist 1,1,1,1,0 -hm 2000,2000,2000,2000,50

echo -e "\nVisualization of ChIP-Seq and RAP-Seq data as histograms"
python HiCPlotter.py -f data/HiC/Mouse/mES.chrX -n mES -r 40000 -chr chrX -o RAP -fh 0 -hist data/HiC/Mouse/GSE46918_pSM33-0hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-1hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-2hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-3hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/GSE46918_pSM33-6hr-Xist_vs_Input.W10000_O7500.bedGraph,data/HiC/Mouse/wgEncodeLicrHistoneEsb4H3k27me3ME0C57bl6StdSig.chrX.bedGraph -hl Xist_0h,Xist_1h,Xist_2h,Xist_3h,Xist_6h,H3K27me3_0h -pi 0 -ptr 0 -fhist 0,1,1,1,1,0 -hmc 4 -sn 0

echo -e "\nVisualization of ChIP-Seq as histograms, ChIA-Pet as arcs and Polycomb domains as tiles"
python HiCPlotter.py -f data/HiC/Mouse/mES.chr3 -n mES -chr chr3 -o Bhlhe22 -r 40000 -s 400 -e 475 -a data/HiC/Mouse/mESC_SMC_ChIPPet.bed -al SMC -hist data/HiC/Mouse/GSM747534_chr3.bedGraph,data/HiC/Mouse/wgEncodeLicrHistoneEsb4H3k27me3ME0C57bl6StdSig.chr3.bedGraph -hl CTCF,H3K27me3 -pi 0 -ptr 0 -t data/HiC/Mouse/mm9_Polycomb_domains.bed -tl Polycomb -tc 00CCFF -ac B4B4B4 -fh 0

echo -e "\nVisualization of 4C data as histograms and Enhancers as tiles with text"
python HiCPlotter.py -f data/HiC/Mouse/mES.chr6 -n mES -chr chr6 -r 40000 -o Digit.vs.GT -s 1295 -e 1338 -hist data/HiC/Mouse/GSM1524258_segToFrag_4C_Digits_WT_E12-5_HoxA13_smoothed_11FragsPerWin.bedGraph,data/HiC/Mouse/GSM1524259_segToFrag_4C_GT_WT_E15-5_HoxA13_smoothed_11FragsPerWin.bedGraph -hl Digits,GT -fhist 1,1 -fh 0 -pi 0 -hm 1500,1500 -pcd 1 -pcdf data/mES_domains_mm9.bed -sn 0.4 -t data/HiC/Mouse/LonfatPrimers.bedGraph -tl Enhancers -tt 1 -ptr 0

echo -e "\nHighlighting selected loci"
python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -chr chr10 -r 25000 -pi 0 -fh 0 -o High -a data/HiC/Human/GM12878.Rad21.bed data/HiC/Human/K562.Rad21.bed -al ChIA-PET Rad21 -s 3000 -e 3500 -pcd 1 -pcdf data/HiC/Human/GM12878_Arrowhead_domainlist.bed data/HiC/Human/K562_Arrowhead_domainlist.bed -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -hl DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq -fhist 0,0,1 0,0,1 -pptd 1 -high 1 -hf data/HiC/Human/highlight.bed 

echo -e "\nAnnotating certain interactions"
python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/KBM7-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HMEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 KBM7 K562 HUVEC IMR90 HMEC NHEK -r 25000 -pi 0 -fh 0 -o Loops -chr chr10 -peak data/HiC/Human/GSE63525_GM12878_replicate_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_KBM7_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_K562_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_HUVEC_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_IMR90_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_HMEC_HiCCUPS_looplist.bed data/HiC/Human/GSE63525_NHEK_HiCCUPS_looplist.bed -s 3600 -e 3675  -ptr 0

echo -e "\nPlotting whole chromosome"
python HiCPlotter.py -f data/HiC/Human/GSM1184323-HiCMYZ-Tcell-Rad21WT-R1.mm9.NA.L-1400000-wDiag-noSS-iced.2.matrix data/HiC/Human/GSM1184321-HiCMYZ-Tcell-Rad21KO-R1.mm9.NA.L-1400000-wDiag-noSS-iced.2.matrix -n Tcell_WT Tcell_Rad21KO -chr Genome -r 1400000 -o Tcell -pi 0 -ptr 0 -wg 1 -hmc 5 -fh 4

echo -e "\nPlotting 5C data"
python HiCPlotter.py -f data/5C/GSM873926_mESCs-female-PGK12.1-day2-Replicate1.txt data/5C/GSM873932_femaleXO-mESCs-DXTX-replicate-1.matrix.txt data/5C/GSM873924_female-MEFs-replicate-1.matrix.txt -n mESC mESC_XO MEF -fh 8 -chr chrX -o 5C -sn 2 -pi 0 -rb 1 -e 300 -hmc 5


echo -e "\nBasic usage"
python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 HUVEC NHEK IMR90 -chr chr10 -r 25000 -s 3000 -e 3500 -o Fig1 -fh 0

echo -e "\nAdding tracks"
python HiCPlotter.py -f data/HiC/Human/GM12878-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/HUVEC-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/NHEK-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/IMR90-chr10_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 HUVEC NHEK IMR90 -chr chr10 -r 25000 -s 3000 -e 3500 -o Fig2 -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseHuvecRawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneHuvecCtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseNhekRawRep2.chr10.bedGraph,data/HiC/Human/wgEncodeBroadHistoneNhekCtcfStdSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqNhekWaveSignalRep1.bedGraph data/HiC/Human/wgEncodeOpenChromDnaseImr90BaseOverlapSignal.chr10.bedGraph,data/HiC/Human/wgEncodeSydhTfbsImr90CtcfbIggrabSig.chr10.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqImr90WaveSignalRep1.bedGraph -fh 0 -fhist 0,0,1 0,0,1 0,0,1 0,0,1 0,0,1 -hl DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq DNAse,CTCF,RepliSeq -hm 400,600,100 400,600,100 400,600,100 400,600,100 400,600,100 -pcd 1 -pcdf data/HiC/Human/GM12878_Arrowhead_domainlist.bed data/HiC/Human/K562_Arrowhead_domainlist.bed data/HiC/Human/HUVEC_Arrowhead_domainlist.bed data/HiC/Human/NHEK_Arrowhead_domainlist.bed data/HiC/Human/IMR90_Arrowhead_domainlist.txt -t data/HiC/Human/GM12878_18_core_K27ac_dense2.bed data/HiC/Human/K562_18_core_K27ac_dense2.bed data/HiC/Human/HUVEC_18_core_K27ac_dense2.bed NHEK_18_core_K27ac_dense2.bed data/HiC/Human/IMR90_18_core_K27ac_dense2.bed -tl ChromHMM ChromHMM ChromHMM ChromHMM ChromHMM -pptd 1 -high 1 -hf data/HiC/Human/fig2.bed 

echo -e "\nCohesin ChIA-PET interactions coincide with early replication sites"
python HiCPlotter.py -f data/HiC/Human/GM12878-chr15_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr15_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -chr chr15 -r 25000 -s 1800 -e 2250 -o Fig3 -hist data/HiC/Human/wgEncodeUwDnaseGm12878RawRep1.chr15.bedGraph,data/HiC/Human/wgEncodeBroadHistoneGm12878CtcfStdSig.chr15.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/wgEncodeUwDnaseK562RawRep1.chr15.bedGraph,data/HiC/Human/wgEncodeBroadHistoneK562CtcfStdSig.chr15.bedGraph,data/HiC/Human/wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -fh 0 -fhist 0,0,1 0,0,1 -hl DNase,CTCF,RepliSeq DNase,CTCF,RepliSeq -hm 400,400,100 400,400,100 -t data/HiC/Human/GM12878_Enhancer.bed,data/HiC/Human/GM12878_Txn.bed,data/HiC/Human/GM12878_Het.bed data/HiC/Human/K562_Enhancer.bed,data/HiC/Human/K562_Txn.bed,data/HiC/Human/K562_Het.bed -tl Enhancer,Transcribed,Heterochromatin Enhancer,Transcribed,Heterochromatin -a data/HiC/Human/GM12878.Rad21.bed data/HiC/Human/K562.Rad21.bed -al RAD21 RAD21 -ptr 0 -high 1 -hf data/HiC/Human/fig3.bed

echo -e "\nA lincRNA locus exhibits active chromatin formation in K562 cells"
python HiCPlotter.py -f data/HiC/Human/GM-chr19_25kb.RAWobserved_KRnormalizedMatrix.txt data/HiC/Human/K562-chr19_25kb.RAWobserved_KRnormalizedMatrix.txt -n GM12878 K562 -r 25000 -chr chr19 -hist data/HiC/Human/GM12878.DNAse.chr19.2.bedGraph,data/HiC/Human/GM12878.RnaSeq.chr19.2.bedGraph,wgEncodeUwRepliSeqGm12878WaveSignalRep1.bedGraph data/HiC/Human/K562.DNAse.chr19.2.bedGraph,data/HiC/Human/K562.RnaSeq.chr19.2.bedGraph,wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph -hl DNAse,RNASeq,RepliSeq DNAse,RNASeq,RepliSeq -t data/HiC/Human/GM12878_TSS+Trx.2.bed data/HiC/Human/K562_TSS+Trx.2.bed -tl ChromHMM ChromHMM -high 1 -hf data/HiC/Human/region.bed -o Figure4 -s 1100 -e 1302 -hm 300,300,100 300,300,100 -fh 0 -fhist 0,0,1 0,0,1 -ptr 0


echo "\n #### test run finished successfully ####"
