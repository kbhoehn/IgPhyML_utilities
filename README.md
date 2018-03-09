# IgPhyML_utilities

This repo contains scripts for simulating sequences under the HLP17 model, and
reconstructing intermediate sequences from IgPhyML output. These scripts were originally just meant for a limited set of analyses, so they're a little clunky to use. Let me know if need help using them or get results that don't make sense, and I'll look into it.

See README.txt for program usage.

Most current public explanation of HLP17 model:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5419485/

You'll need the most recent version of IgPhyML installed: 
https://github.com/kbhoehn/IgPhyML

You'll also need the Perl programming language and the following Perl 
dependencies installed:

PDL

PDL::LinearAlgebra::Trans


I recommend using cpanm to install them.

If you use these scripts (or part of them) for a publication, please cite:

Hoehn, K. B., Lunter, G., & Pybus, O. G. (2017). A phylogenetic codon substitution 
model for antibody lineages. Genetics, 206(1), 417-427. 
