CRISPRandME
===============================


a proof-of-concept package of optimal gRNA search tool with user-defined reference genomes

Quick-Start
--------
* Clone git repo
`git clone https://github.com/BMES550-CRISPR/CRISPRandME.git`

* Install required python modules
`conda install --file requirement.conda`

* Run test case
` python bwt.py -i ../raw_data/chr22_10k.fasta.gz -q TGCGGAGAGTGTGCGGCTCCAGG -g "['HG01390', 'HG01889', 'HG02455', 'HG02943','HG02982', 'HG03084', 'HG03291', 'HG03514', 'HG03538', 'HG03687', 'HG03950','NA11918', 'NA19209', 'NA19213', 'NA19320', 'NA19380', 'NA19454', 'NA19455','NA19717']"`

`-i: reference genome sequence. 10k bp of chr22 is used for demonstration.`

`-q: gRNA sequence. 20 bp protospace + 3 bp PAM (NGG) for the search of S. pyogenes Cas9 gRNAs.`

`-g: group of individuals in list selected from 1000 Genome Project phage 3 (2,504 individuals involved).`
