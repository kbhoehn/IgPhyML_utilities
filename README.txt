# IgPhyML_utilities

A phylogenetic codon substitution model for antibody lineages
Supplemental code used in ancestral state reconstruction
Kenneth B Hoehn
kenneth.hoehn@yale.edu
9/March/2018


Introduction

This repo contains scripts for simulating sequences under the HLP17 model, and
reconstructing intermediate sequences. These scripts were originally just meant 
for a single publication, so they're a little clunky to use and may fail under 
certain weird inputs. Let me know if you get any odd results and I'll look into it.

The supplied script ancReconstructHLP17.pl is a perl script that calculates the 
most likely marginal ancestral sequence at each node and site, given a number 
of parameters taken from an IgPhyML analysis. It works by reading in parameters 
and files from a config file. See the files in the examples/ directory for 
examples.

simulateHLP17.pl is used for simulations. It's description is in the last 
section.

hotSpotUtils and hotSpotUtils.pm contain the subroutines used by these scripts.
You'll need to set the environment variable PERL5LIB to the directory where
hotSpotUtils.pm and hotSpotUtilsPDL.pm are located (i.e. this directory). To do this
you'll probably need add 

export PERL5LIB=$PERL5LIB:<PATH TO THIS DIR>

to your ~/.bashrc file

For explanation of marginal ancestral reconstruction see Pupko et al (2000).

For more detail into the specific caluclations used here see Boussau and Gouy
(2006) and Boussau et al (2008).

Most current public explanation of HLP17 model:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5419485/

You'll need the most recent version of IgPhyML installed: 
https://github.com/kbhoehn/IgPhyML

You'll also need the Perl programming language and the following Perl 
dependencies installed:
PDL
PDL::LinearAlgebra::Trans

I recommend using cpanm to install them.

If you publish any analyses using these scripts (or part of them), please cite:

Hoehn, K. B., Lunter, G., & Pybus, O. G. (2017). A phylogenetic codon substitution 
model for antibody lineages. Genetics, 206(1), 417-427. 


1. When you run IgPhyML on a data set
#####################################

Use the --ambigfile <ambigfile> flag when you use the HLP17 model. This will 
output a file which specifies the assignments IgPhyML makes on ambiguous 
nucleotides (e.g. N's). Record the name of this file in the config file. For 
the example, run:

igphyml -i example/CH103.20.phy -m HLP17 --motifs WRC_2:0,GYW_0:1 --hotness e,e 
--root V4-59 --run_id hlp --ambigfile example/CH103.ambig

Which will fit an asymmetric WRC/GYW motif model. If you want to try a 
different motif model, it's pretty easy - check out the manual.


2. Complete config file
#####################################

Within the config file you need to specify IgPhyML stats and (re-rooted) tree 
files, location of IgPhyML installation directory, ambiguous character files 
(see point 1.), output directory, output stem, germline root ID, sequence file 
(fasta format), length of the sequences, and whether the tree is rooted (for 
now it should always be). Sequence length is in codons.

The example config file CH103.example.config.txt contains the parameters:

length	81
rooted	1
outdir example
seqfile	example/CH103.20.fa
rootid	V4-59
igphyml	../../../IgPhyML_development/IgPhyML_dev
stats example/CH103.20.phy_igphyml_stats.txt_hlp
tree example/CH103.20.phy_igphyml_tree.txt_hlp.figtree.reroot
ambigfile example/CH103.ambig 
stem	CH103.example

To run this you need to replace the current value of "igphyml" with the 
directory on your computer. If you've followed the installation instructions in 
the IgPhyML manual, you can find this easily by entering the command 'which 
igphyml' and omitting the 'src' at the end. 


3. Run script
#####################################

The basic operation is

perl ancReconstructHLP17.pl <config file> <overwritten parameters>

For example:

perl ancReconstructHLP17.pl CH103.example.config.txt


To do these analyses on your own files, you'll need to specify these parameters 
according to your own analysis. You can also overwrite any of these parameters 
by adding a - to them and specifying their value when you run the program. For 
instance, to use a different output stem, run:

perl ancReconstructHLP17.pl CH103.example.config.txt -stem CH103.differentfile


If all goes well a lot of text will print out showing the settings of the 
reconstruction run. Check to make sure they are what you want. After the 
initial deluge of information a slower trickle will follow, showing first the 
lower partial likelihoods of each subtree as they are calculated. Next will 
follow the upper partial likelihoods of each upper tree, specified by their 
subtaxa. Next will come information about the marginal ASRs of each node. Part 
of sanity checking this caluclation is printing out the likelihood at each 
node, which should be the same as that calculated in the IgPhyML stats file, 
given the latter is rounded to two decimal places. Finally you'll get the 
likelihood calculation at the root node compared to that from the stats file. 
If this, or any of the previously mentioned likelihood values spat out during 
the ASR stage (but not lower/upper partial likelihoods) are different from that 
in the stats file, something is wrong. Email me if you can't figure it out.

One thing to keep in mind is that the reconstructed codons are the most likely 
codon reconstructions at each node (marginal), rather than the most likely sequence 
of mutations across the tree (joint), so be careful comparing the reconstructions 
between nodes. Differences in codons between nodes don't necessarily indicate a 
substitution likely occurred between them. I'd also recommend directing the output 
into a log file, because part of the terminal output describes confidence in each 
site's assignment. Email me if you need an explanation of that.

For 20 taxa, as in the example file, this doesn't take long, but if you have a 
very large dataset you may need to watch these results trickle in for a while.

Finally, this will result in two output files:

<outdir>/<stem>.MLcodons.fa

and 

<outdir>/<stem>.MLaas.fa

These have the most likely marginal reconstructions at each node for codons and 
amino acids, respectively. Their ids contain:

<Cladistic degree>;<Subtaxa1,..,SubtaxaN>;<Divergence from root>

Using these peices of information you can figure out where on the tree these 
sequences fit.


Simulations
#####################################

Simulations use a similar configuration file as reconstructions. See sim.config for an example.

perl simulateHLP17.pl <config file> <overwritten parameters>

for example (once you properly set the 'igphyml' flag - see section 2):

perl simulateHLP17.pl sim.config

This will output 10 simulated lineages using the parameters specified in 
example/CH103.20.phy_igphyml_stats.txt_hlp and the tree in 
example/CH103.20.phy_igphyml_tree.txt_hlp.figtree.reroot. 

It is possible to do fully context-dependent simulations under the HLP17 model, but these
are a bit more complicated. You can get an example in the code from File S2 from: 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5419485/

References
#####################################

Boussau B, Blanquart S, Necsulea A, Lartillot N, Gouy M. 2008. Parallel 
adaptations to high temperatures in the Archaean eon. Nature 456:942–945.

Boussau B, Gouy M. 2006. Efficient Likelihood Computations with Nonreversible 
Models of Evolution. Syst. Biol. 55:756–768.

Pupko T, Pe I, Shamir R, Graur D. 2000. A fast algorithm for joint 
reconstruction of ancestral amino acid sequences. Mol. Biol. Evol. 17:890–896.