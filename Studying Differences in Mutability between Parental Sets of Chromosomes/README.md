# Studying-differences-in-mutability-between-parental-sets-of-chromosomes
Different genomes (i.e. sets of chromosomes) differ at many positions in sequence and genomic structure. Both characteristics, especially the latter (chromatin state), are known to impact mutation rates in our tissues, but how much of a difference does it make between parental sets of chromosomes?

### Planning
__Step 1__: PCA of SNPs from the 1000 Genome project (SNPs on exons that are present in PCAWG( is the buca patients etc …)) and projections of SNPs from PCAWG onto the principal components calculated on the 1000 Genome project. -> Allows identification of mixed ancestry samples in PCAWG. <br>
__Step 2__: For 1 mixed ancestry sample, phase the SNPs to an ancestry by blocks of 500,000 bases (you may have to play with the window size later). <br>
__Step 3__: For this sample, phase the SNVs to the SNPs and thus to the ancestry. <br>
__Step 4__: Generate graphs or statistics on the somatic mutation rate on each parental copy for this sample. <br>
