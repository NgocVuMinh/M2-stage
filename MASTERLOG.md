# MASTER LOGFILE

## Mon 19 Jan
Preliminary resources:
- [Neurology course by Ninjia Nerd](https://youtube.com/playlist?list=PLTF9h-T1TcJgx3OFachdjHPMX6VE4VDS1&si=exX_IpoPaBiapoTe)
- [Neurogenesis in the mammalian brain by Susannah Hannaford](https://www.youtube.com/watch?v=B2RINOAeONw)
- [Neuroscience: Exploring the brain, 3rd Edition, Chapter 23. Bear, Connors & Paradiso](https://seti.net/Neuron%20Lab/NeuronReferences/Neuroscience%20-%20Bear.pdf#page=727). Video lecture: [OSSM Neuro Chapter 23 - Neurogenesis by Brent Richards](https://www.youtube.com/watch?v=pSb9iHVnj_8)

Baffet lab's papers:
- [Coquand 2024, NatCell](papers/Coquand2024_NatCell.pdf) and [PR report](papers/Coquand2024_NatCell_PR.pdf). 
- [Lindenhofer 2024, NatCell](papers/Lindenhofer2024_NatCell.pdf) and [PR report](papers/Lindenhofer2024_NatCell_PR.pdf)


## Thu 22 Jan - Fri 23 Jan
Meeting with OS: Overall: 
- Phase I: Exploring / analyzing the HDBR dataset by Lindsay 2016 of normal developing brain samples --> benchmark / reference transcriptomic signatures for each step of the brain development (early embryonic?) of each brain compartment (sounds ambitious)
- Phase 2: Analyzing Ludovica's celebellar organoids ie mapping to the phase I dataset
So for the second step to be good the first step should be good, hope so.

Data is ready, QC and alignment by Loic, will need to ask for QC reports and alignment reports to see them myself.

2 streams of analyses are DEGs and AS. Loic will handle AS mostly? Will confirm later with OS.

Phase I:
- Analyze at the highest level first: mid, hind, fore --> stratify further if needed
- 2 QC steps from the perspective of the dataset: (1) reproduce the PCA in the HDBR paper (2) deconvoluting existing samples:
    (1) PCA: maybe 3 different PCAs: with and without the unidentified samples, and unident samples alone (to see if they form any clusters on their own). It looks like they could serve as a blind testing set if i come up with a way to deconvolute then predict sample type / dev stage based on transcriptional profile (unlikely but we'll see).
        The PCA results of cleaned samples will be used to re-stratify the samples from separate time points to grouped dev stage (eg early=4-9PCW) but be mindful that the grouping might (very likely) differ across brain compartments --> consider doing separate PCA on individual brain comps first, then aggregates as the variance will start small, thus easier to identify trajectory in PCA at first. ideally, in all-comp PCA, PC1-3 should explain both dev stage and brain comp, and single-comp PCA should exlain dev stage.
        The clusters identified in PCA after careful calibrations (number of top variable genes or things like that + literature justifications) will be used for 2 things: (1) grouping time points then (2) DEGs (expect lots of back and forth / iterative pairwise comparisons).

    (2) Deconvolution of samples: since all are bulk seq, need to know whether they are pure / heterogenous. OS mentioned something like estimating the % of cells of each cell type in a sample?

Presentation next tuesday 27 Jan.

Loic is not done yet with the raw data (need to merge them?), he said he needs to try 3 different protocols to see which one produced the highest % of mapped reads (if I understood him correctly), then picked that one. So I'll probably not have any interesting things to show on my 1st presentation next tuesday 27 Jan. 

In the meantime i work on the metadata. Talked to OS, he's okay: I need to compare the metadata they published (supp tab1 found on HDBR website and the one they actually deposited with their data). I got the latter from Loic already but havent looked at it. Looked into tab1, found it's not consistent entirely with the paper's table 1, so will need to check that again (use `scr/metadata.R`). Will compare with Loic metadata.


REF to follow: https://www.sciencedirect.com/science/article/pii/S2211124721008597
https://academic.oup.com/nargab/article/2/3/lqaa078/5909519


Also got my soroban access. Will need to re-learn the thing on Monday (with Loic?). Got a tutorial docs from Lucas too (in `docs/tut_soroban.docx`?).


## Thurs 29 Jan

Should have written this more consistently.

Before the meeting I did metadata analysis. 

1st meeting was good. People mostly asked questions about the sample preparation. I didnt know abortion was that complicated. 

I did batch correction first before 

Justifications to remove samples that come from batches with lower than 5 samples: https://academic.oup.com/biostatistics/article/17/1/29/1744261#137601345.

--> remove and try ComBat-seq





