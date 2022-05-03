#!/usr/bin/env python
# coding: utf-8

# # 4. Immunogenicity analysis
# 
# We uses the method of removing and/or reducing potential T-cell epitopes, as an
# approach to the management of the immunogenicity of biologics. The protein sequence is
# scanned in silico, for sequences that have a strong binding signature for a family of 50 MHC
# Class II receptors , whose alleles cover 96 – 98% of the human population. The presented4
# histograms for each variable region sequence, show the average (for the n positively-testing
# MHC II alleles) of epitope strength at each position as a percentage for all epitopes above a
# threshold of 20%. At each position in the sequence, the number of alleles scoring above the
# threshold is shown above the histogram at that position. The epitopes of most concern for the
# antibody’s immunogenicity are therefore those that have not just the highest average score per
# allele (as shown by the histogram), but which also score above the threshold across more
# alleles, since these epitopes are more likely to engender an immune response in a larger
# fraction of the patient population.
# 
# Experience using in silico algorithms of this kind in conjunction with laboratory immunogenicity
# assays has shown that epitopes below this threshold do not generally contribute significantly to
# the protein’s immunogenicity. The number of alleles, the affected alleles and their individual
# scores are also listed in the detailed analyses below each histogram figure.
# 
# The raw immunogenicity score quoted is the total over all epitopes above the threshold for all
# affected alleles. The normalized immunogenicity score is this raw score divided by the
# sequence length, and represents epitope strength per unit sequence to enable comparisons of
# protein sequences of different lengths.
# 
# The absolute magnitudes of these scores are somewhat arbitrary, but they have value as
# comparative metrics. It has been shown that human serum proteins generally display an
# immunogenicity potential that is inversely proportional to their abundance in serum . Proteins5
# that are found at very low concentrations in serum, like erythropoietin, can have normalized
# scores above 80%. By contrast, very abundant human serum proteins like albumin and
# immunoglobulins typically have normalized scores in the 35 - 50% range.
# 
# 
# ## immunogenicity analysis - heavy chain
