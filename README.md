# EPhIMM
**E**xpress **Ph**ylogenetic **I**nference based on **M**ultiple **M**arkers

**EPhIMM** is a shell script for express phylogenetic analysis of prokaryotic genomes.


## Dependencies

* **esl-sfetch** from Easel for <https://github.com/EddyRivasLab/easel>
* **mafft** <https://mafft.cbrc.jp/alignment/software/> for multiple sequence alignment
* **fasttree** <http://www.microbesonline.org/fasttree/> for tree inference
* **Java 8** for autoedit of alignment (typically installed in Linux ditributions)

Install some of them with `conda install mafft fasttree easel`

If You start from nucleotide genomic sequences You will need **prodigal** <https://github.com/hyattpd/Prodigal>. It could be installed with *conda* or Your linux distribution's package manager like *apt*, *yum* etc.


## How it works

Extract HMM profiles for marker genes from `markers.tar.gz` using `tar xzf markers.tar.gz`.

Run EPhIMM.sh in directory containing files with protein sequences. E.g. if You downloaded EPhiMM to `/home/user/software/EPhiMM/` and data located in `/home/user/data/mygenomes`:

`user@pc:~/ cd /home/user/data/mygenomes`
`user@pc:/home/user/data/mygenomes$ /home/user/software/EPhiMM/EPhiMM.sh`

Default output folder is `output-DATE`, where DATE is current date. You can provide different output folder using flag **-o**, e.g. `EPhiMM.sh -o MyOutput`.

If You work with nucleotide genomic sequences use flag **-p** to predict CDS with Prodigal.

After marker genes search. fetching and alignment You will be asked for manual or automatic alignment editing. In case of automatic columns containing more than 50% gaps will be removed from the alignment using java tool *aleditor.jar*. In case of manual editing open file `all.aligned.faa` in the output folder, make changes and save it. You will be prompted for relative path to edited file, e.g. if You saved file as edited.version.faa in final directory type `final/edited.version.faa`.

You can check `genomes.stats.txt` and `hmm.stats.txt` for missing and duplicated genes. Resulting alignment and tree will be stored in `final/` directory.


## Markers

Marker genes should be provided as HMM profiles.

Initial set containing 41 single copy marker genes was taken from CheckM software (Parks et al., 2015), 2 markers from the original set(TIGR00344 and TIGR00422) were excluded due to multiple hits.

If You like to use Your own set of markers put `.hmm` files in a directory, e.g. /home/user/MyMarkers, and run the tool providing path to the markers with flag **-m**: `EPhIMM.sh -m /home/user/MyMarkers/`.

Default e-value for hmmsearch is 1e-6, You can set it with glag **-e**: e.g. `EPhIMM.sh -e 1e-10`.


## References

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490. 
