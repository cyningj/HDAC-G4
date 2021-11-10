# RNA-seq data analysis of MDA-MB-231 cells treated with HDAC/G4 dual-targeting compound

1.G4Hunter prediction

We extracted the sequences of 2000 bp upstream of TSS of the Differentially expressed genes (DEGs) for G4 prediction. First, we used the annotation file 'gene_encodev36' to get DEGs' TSS, gene ID, sense/antisense strand and chromosome information. Second, we used getfasta (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html) to extract the sequences of 2000 bp upstream of TSS. Then, we predicted the G4H of all the sequences using a sliding window approach with a window length of 60 nt and a step of 5 nt. by G4Hunter ([http://bioinformatics.ibp.cz](http://bioinformatics.ibp.cz/).).



2.G4 annotation

The successive windows with G4H higher than 1.2 were grouped and returned as pG4 regions (pG4r) containing one or more pG4. Then two sets of detected pG4r were merged in order to detect and remove redundancy, and G4 score was counted according to the number of optimized pG4r using the script 'G4Detect-j.py' which we modified based on the script 'G4Annotation.py' from [G4HumanTranscriptome](https://github.com/UdeS-CoBIUS/G4HumanTranscriptome). 100 genes were randomly selected from the human genome and were run the same protocol for the G4 background score. DEGs scored over the background in G4 score were defined as G4-related DEGs.
