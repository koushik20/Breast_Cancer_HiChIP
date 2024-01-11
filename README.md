# Breast_Cancer_HiChIP
HiChIP to identify risk genes associated with breast cancer susceptibility, shedding light on the complex regulatory interactions underlying the disease.
Integrating HiChIP data of 5 breast cell lines (MCF7, MCF10A, T47D, MDAMB231, Htert-HMEC) and analyzing H3K27ac-marked enhancer-gene interactions, we can uncover non-coding genetic variants and target gene regulatory interactions involved in breast cancer risk loci.

Step-1: HiC-Pro to analyze raw fastq reads and generate valid read pairs
Step-2: FitHiChIP statistical analysis by integrating H3k27ac peaks to find significant cis interactions
Step-3: diffloop to find differentially expressed loops between the cell lines
Step-4: Integration of various datatypes such as TWAS, eQTL, credible causal variants etc
