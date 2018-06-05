# kmt2d

This work relates to my thesis project. It is almost entirely R-based, with a few shell scripts.

KMT2D is a histone methyltransferase. It was the 13th most mutated gene in a pan-cancer analysis and mutation rates of 89% in follicular lymphoma have been reported. The hypothesis I was investigating was whether it has a role in DNA repair. To decipher whether this was the case, I had the following data:

- 600 lymphoma samples with Strelka variant calls in variant call format (VCF) and mutation annotation format (MAF) files
- Obtain list of filepaths in txt file, from which to pull out mutational calls of user-defined transcript
- Aligning of snv and indel call files for each sample and definition of all samples as KMT2D wildtype or mutant

and took the following approach:

## 1) Overview of data characteristics 
a) St Jude Pecan Cloud gene mutation visualisation tool 
  i) Script to wrangle data into format accepted by tool
  ii) Compare to mutations listed in COSMIC and those from St Jude paediatric cancer samples
b) Stacked bar graph of number of KMT2D mutated samples as a proportion of all samples
  
  
## 2) Visualisation of KMT2D mutational status correlation with SNV/indel counts
a) Correlation of proportion of samples with mutated KMT2D to median somatic alterations by lymphoma type
b) Distribution densities
c) Beeswarm plots
d) Boxplots
e) Provean online tool used to analyse which mutations are likely to be deleterious and which not (this line of questioning was abandoned when I discovered that synonymous mutations are often cancer driver genes (REF) for various reasons such as codon bias, and that Provean delineated all synonymous mutations as wildtype)



## 3) Comutation or synthetic lethality with other DNA repair genes
a) PARP1
b) BRCA1
c) Maftools summaries: concatenated Maf files using an R script to print all file names within a lymphoma type and a shell script


## 4) Mutational signature analysis: is there a known signature associated with KMT2D mutational status, and if not, is there a novel one?
a) Mutational trinucleotide context extraction
b) Distillation into 30 COSMIC-defined mutational signatures (some with known aetiology)
c) Heatmaps for both drawn with KMT2D status aligned along column labels






### Problems I ran into & how I solved them
- Consistency between samples and data available: Maf or Vcf or both file types available
- FFPE samples
- Obtaining filepaths to import into RStudio owing to inconsistent file naming 
- For each sample, the indels and the snv calls were stored in separate files. To get the mutational status of an entire sample, it was necessary to find unique identifier strings within the filenames that would be shared by snv and indel file and no others
- When loading Vcf files into RStudio, depending on the sequencing run and settings, differing numbers of columns were present, causing alignment and rbinding issues
- Calls in Vcf files in which the transcript sought was not the *main* transcript listed, but rather in a column containing a mix of all transcripts that could be affected and how. This required extraction of the clause detailing the effects to that transcript and relocation to the relevant columns.
