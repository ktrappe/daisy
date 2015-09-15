# Daisy
#### Horizontal Gene Transfer Detection by Mapping Sequencing Reads

Daisy is a pipeline for horizontal gene transfer (HGT) detection from
sequencing data. It requires sequencing data from the HGT organism and
reference sequences from the acceptor/recepient genome and the donor genome.

Daisy is a pipeline written in Python that uses Samtools and Bamtools for SAM
file processing and extraction of unmapped reads, the C++ SeqAn tools Yara and
Gustaf for mapping, and contains a Python based evaluation routine.

Contact: kathrin.trappe@fu-berlin.de

## Installing Daisy
Daisy needs to have Python installed with a version supporting `subprocess32`.
The easiest way to get Daisy is to
download **daisy.py** and the script **hgt_eval.py** and place both scripts in your
`~/bin/` directory.
Daisy depends on the following established open-source tools which have to
be installed either in your `~/bin/` as well or be globally avaible for all
users on your server.

###### Samtools, Bedtools
To install Samtools and Bedtools, please follow the installation guides given
on
http://samtools.sourceforge.net/
and
http://bedtools.readthedocs.org/en/latest/


###### Yara, Stellar, Gustaf, SAK
Precompiled binaries (Linux 64-bit, Windows, Mac OS X) of Yara, Stellar, Gustaf,
and SAK can be downloaded from the SeqAn projects download page:
http://www.seqan.de/downloads/projects.html

All tools are distributed with SeqAn - The C++ Sequence Analysis Library (see
http://www.seqan.de). To build them yourself, you can check out the latest
developer version of SeqAn on:
http://github.com/seqan/seqan/tree/develop/
Follow the installation guides given on
http://seqan.readthedocs.org/en/master/Tutorial/GettingStarted/LinuxMakefiles.html


###### Phage database
As an additional step, Daisy maps the reads against a phage database, and flags HGT
candidates having relevant hits. We recommend using the phage database available from
http://www.ebi.ac.uk/genomes/phage.html

## Using Daisy
##### Example
Download the folder "data/example" and make sure you have all tools ready.
The example run is a subsample from the simulated data set below. The reads
stem only from the inserted phage sequence plus 2000bp surrounding sequence.
The run takes only a few minutes. Within the example folder, you can run it via
```
python ~/bin/daisy.py -r1 Ecoli_K12_mod_HPylori_1322000-1350000_mod_1115289-1147285.1.fa
                      -r2 Ecoli_K12_mod_HPylori_1322000-1350000_mod_1115289-1147285.2.fa
                      -ar ../Ecoli_K12.fa -dr ../Helicobacter_pylori_ML1.fasta
                      -a "gi|170079663|ref|NC_010473.1|" -d "gi|766541424|dbj|AP014710.1|"
```
If you have downloaded the phage database, add its path to the program call
using the '--phage_ref' option (e.g. --phage_ref phage_all.fasta).
The produced result files should be the same as the corresponding *gold* files
provided in the folder.
##### H. pylori data set
The H. pylori data set is the complete simulated data set evaluated in the paper.
From within the HPylori folder, you can re-run it via
```
python ~/bin/daisy.py -r1 Ecoli_K12_mod_HPylori_1322000-1350000_mod.1.fasta
                      -r2 Ecoli_K12_mod_HPylori_1322000-1350000_mod.2.fasta
                      -ar ../Ecoli_K12.fa -dr ../Helicobacter_pylori_ML1.fasta
                      -a "gi|170079663|ref|NC_010473.1|" -d "gi|766541424|dbj|AP014710.1|"
                      -new
```
Daisy checks for the presence of already computed files. So if you omit the `-new` parameter, Daisy will recognize the existing files and run through without changing results.
Use the `-task` parameter to assign job names. You can also specify each pipeline step to be run or not run separately (see help message).

## Output Formats
Daisy currently supports the VCF output format for reporting HGT candidates
meeting the pre-defined threshold.
Additionally, all HGT candidates together with their sampling results are
written to a TSV file.

##### Variant Call Format (VCF)
The output is according to VCF 4.2. We report the single HGT boundaries as
inter-chromosomal translocations as SV type, connect the boundary pairs via
identical IDs and introduce the event tag 'HGT'.

See http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41
for information about the VCF file format specifications.
