## From Raw Data to the OTU Table 
This wiki is based on [Frédéric Mahé's metabarcoding pipeline](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline), which I think is a great and simple pipeline for building the OTU table. Here, I adapted the commands to my needs.

## Requirements
For this pipeline we need [vsearch](https://github.com/torognes/vsearch), [cutadapt](https://github.com/marcelm/cutadapt/) and [swarm](https://github.com/torognes/swarm). The final python script for building the OTU contingency table only works with `python 2.7`

Here are the versions and authors of the used programs (as of May 2019):

    ## vsearch v2.10.3_linux_x86_64, Rognes et al. 2016
	## cutadapt version 1.18, Martin 2011
	## Swarm 2.2.2, Mahe et al. 2015

## Merge Raw data
The Cologne Center for Genomics sequenced our samples on a MiSeq platform with 250bp paired reads. As we amplified the ITS1 region of the oomycetes, which is  ca. 250-350bp long, we expect quite a long overlap between the paired reads.

We merge our reads with `vsearch`, like this:

```sh
VSEARCH=$(which vsearch)
THREADS=46
ENCODING=33
FORWARD="00_KD09-Oomycota-150418_S1_L001_R1_001.fastq.gz"
REVERSE="00_KD09-Oomycota-150418_S1_L001_R2_001.fastq.gz"
OUTPUT="01_Oomycota.merged.fastq"

"${VSEARCH}" \
    --threads ${THREADS} \
    --fastq_mergepairs ${FORWARD} \
    --reverse ${REVERSE} \
    --fastq_ascii ${ENCODING} \
    --fastqout ${OUTPUT} \
    --fastq_allowmergestagger \
    --quiet 2>> ${OUTPUT/.fastq/.log}
```
Looking at the generated log file, we see some statistics after the merging (this output is from the Spring-Autumn_2018 samples, which is not published yet):

      14888229  Pairs
      13223430  Merged (88.8%)
       1664799  Not merged (11.2%)

    Pairs that failed merging due to various reasons:
      105145  too few kmers found on same diagonal
      264134  multiple potential alignments
      576099  too many differences
      718516  alignment score too low, or score drop to high
         905  overlap too short

    Statistics of all reads:
      238.56  Mean read length

    Statistics of merged reads:
      303.59  Mean fragment length
       78.47  Standard deviation of fragment length
        0.26  Mean expected error in forward sequences
        0.43  Mean expected error in reverse sequences
        0.15  Mean expected error in merged sequences
        0.21  Mean observed errors in merged region of forward sequences
        0.41  Mean observed errors in merged region of reverse sequences
        0.62  Mean observed errors in merged region

Nearly 90% of all reads were merged, and the contigs have the expected length of ca 300bp +- 70bp, that's good.

## Demultiplexing, primer clipping, sample dereplication and quality extraction

We have 81 samples, all of them tagged with an individual barcode. These tags can be used to demultiplex our run, so we have seperate files for each sample. The file `Oomycota_Tags.tsv` looks something like this:

     K232_Bl_Fr	GCTTCTAG	GATACTGC
     K232_Tot	GCTTCTAG	GCTGTGAT
     K439_Bl_Fr	GCTTCTAG	TGGAGCTT
     K439_Tot	GCTTCTAG	TTACGCGA
     K232_Bor	TAGCTCTG	TGACAGGT
     K232_Det	TAGCTCTG	AGGCTGTA
     K439_Bor	TAGCTCTG	TCGCTGAA

The first column gives the name/ID of the sample, the other two the forward and reverse barcode, respectively.

Then we use `cutadapt` to demultiplex the samples. `cutadapt` first searches for contigs containing the forward tag, then the forward primer, then the reverse tag, and finally the reverse primer. This way it doesn't matter that some samples have the same e.g. forward tag, because `cutadapt` takes the combination of both tags into account. 

Because we are dealing with merged contigs ("single end" so to say), we need the **reverse complement** of our **reverse primer and tags**. Nevertheless, to avoid losing data, we also search within the **reverse complemented reads**. Quite complicated, but still better than accidentally losing data :) I will go into more detail in the script. So, we use following primer sequences:

    Forward: GCGGAAGGATCATTACCAC
	Reverse: GCTCGCACAHCGATGAAGA <- This is already reverse complemented
	
With that, we can start the demultiplexing, primer clipping and quality extraction:

```sh
INPUT="01_Oomycota.merged.fastq"
TAGS="Oomycota_Tags_ReverseRevComp.tsv"
PRIMER_F="GCGGAAGGATCATTACCAC"
PRIMER_R="GCTCGCACAHCGATGAAGA"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F}))
MIN_R=$(( ${#PRIMER_R}))

# Define binaries, temporary files and output files
CUTADAPT="$(which cutadapt) -j 46 --discard-untrimmed -e 0 --minimum-length ${MIN_LENGTH}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTQ=$(mktemp)
TMP_FASTQ2=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT=$(mktemp)
QUALITY_FILE="${INPUT/.fastq/.qual}"
  
# Reverse complement fastq file
"${VSEARCH}" --quiet --threads 46 \
             --fastx_revcomp "${INPUT}" \
             --fastqout "${INPUT_REVCOMP}"

while read TAG_NAME TAG_SEQ RTAG_SEQ; do
    LOG="${TAG_NAME}.log"
    FINAL_FASTA="${TAG_NAME}.fas"

    # Trim tags, forward & reverse primers (search normal and antisens)
    cat "${INPUT}" "${INPUT_REVCOMP}" | \
        ${CUTADAPT} -g "${TAG_SEQ}" -O "${#TAG_SEQ}" - 2> "${LOG}" | \
        ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2>> "${LOG}" | \
        ${CUTADAPT} -a "${RTAG_SEQ}" -O "${#RTAG_SEQ}" - 2>> "${LOG}" | \
        ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_F}" - 2>> "${LOG}" > "${TMP_FASTQ}"

    # Discard sequences containing Ns, add expected error rates
    "${VSEARCH}" \
        --quiet --threads 46 \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --relabel_sha1 \
        --eeout \
        --fastqout "${TMP_FASTQ2}" 2>> "${LOG}"

    # Discard sequences containing Ns, convert to fasta
    "${VSEARCH}" \
        --quiet --threads 46 \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --fastaout "${TMP_FASTA}" 2>> "${LOG}"
    
    # Dereplicate at the study level
    "${VSEARCH}" \
        --quiet --threads 46 \
        --derep_fulllength "${TMP_FASTA}" \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${FINAL_FASTA}" 2>> "${LOG}"

    # Discard quality lines, extract hash, expected error rates and read length
    sed 'n;n;N;d' "${TMP_FASTQ2}" | \
        awk 'BEGIN {FS = "[;=]"}
             {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
        tr -d "@" >> "${OUTPUT}"
done < "${TAGS}"

# Produce the final quality file
sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
    uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean
rm -f "${INPUT_REVCOMP}" "${TMP_FASTQ}" "${TMP_FASTA}" "${TMP_FASTQ2}" "${OUTPUT}"
```
This might take a while, even with 46 threads. Grab a coffee and have a break in the meantime :)

The script produces 81 fasta files, which look like this (example from K129_Bod.fas):

     >8563321f7bb94c43c2f46697b3923d1c5324b2e4;size=1081
     ACCTAAAAACTTTCCACGTGAACCGTTGAACAAAGATTTTGTGCTTGTCTAGGCTTTTAGTTTGGACTTGCTAATCGAAGGTAGTCGTAAGATTACCGATGTTCTTTCAAACCATTTAATATTTACTGATCCTTACTCTGAAGACGAAAGTCTACAGTTTTAATCCACAACAACTTTCAGCAGTGGATGTCTAG
     >0fbb3376b274ce5299d823e11014741e831677ab;size=796
     ACCAAAAAAAACTTTCCACGTGAACCGTTGTAATTATGTTCTTGTGCTTTCCTTCGGGAAGGCTGAACGAAGGTGAGCCGCTTTATTTTTTGTATCGTGGCTTGCCGATGTACTTTTCAAACCCATTTACTTAATACTGAACTATACTCCGAAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAG
     >31d85ee4271b9553524f8fbf23e0f43e55e9919c;size=783
     ACTCCAAAAAACTTTTTCTAGTCCATATCTTAAATGATAGTCTGATTTGAATTGCAGCCCTGGTGGCGGATTCATTTTTAAAATCACAATGACCGTAGACGATGGATGACTTG

The reads are assigned a unique identifier (some cryptic ID with lots of numbers and letters), and identical reads are collapsed and counted, represented by the `size=` field in the ID. This means that e.g. the first read - or rather this particular sequence - in the above example occurs 1081 times within the sample K129_Bod.

For every .fasta file, there will be a corresponding .log file. It first gives you the command line parameters, then a summary of the processed reads. This .log file is "additive", so the first summary in each file gives you the number of reads containing the **forward tag**:

```
=== Summary ===

Total reads processed:              26,446,860
Reads with adapters:                   635,938 (2.4%)
Reads that were too short:                 101 (0.0%)
Reads written (passing filters):       635,855 (2.4%)

Total basepairs processed: 8,029,113,520 bp
Total written (filtered):    176,574,010 bp (2.2%)
```

This is of course a small number compared to the total amount of reads processed. In this example, these 635.938 reads will then be passed to the next summary, which gives you the reads containing the **forward primer**:

```
=== Summary ===

Total reads processed:                 635,855
Reads with adapters:                   426,533 (67.1%)
Reads that were too short:                  72 (0.0%)
Reads written (passing filters):       426,461 (67.1%)

Total basepairs processed:   176,574,010 bp
Total written (filtered):    110,138,788 bp (62.4%)
```

Then these 426.461 reads will be passed, and so on.

The last two summaries will give you a warning which looks something like this:

```
Bases preceding removed adapters:
  A: 99.7%
  C: 0.0%
  G: 0.2%
  T: 0.1%
  none/other: 0.0%
WARNING:
    The adapter is preceded by "A" extremely often.
    The provided adapter sequence could be incomplete at its 3' end.
```

No need to worry about that, this phenomenon makes perfect sense in our case. 
Why? Well, we designed the primers on conserved regions of the ITS1 sequence, so of course in most cases there will be a particular base preceding the primer. Same for the reverse tags: they are attached to the primer, which has a fixed sequence and a fixed last base preceding the tag. 

We also get a quality file (`01_Oomycota.merged.qual`) listing the expected error rates for all unique sequences from all samples together as well as the sequence length in increasing order:

     01ac17146b57dec5572d1c97f80259d65644146a	0.0025	32
     01f671b236010bc9bab4b5f91d2809bfce27b408	0.0025	32
     027bb0054262fcc9f1ef8718ad740bf61e761392	0.0025	32
     033b189fd804271ef70f338bc41fb7ddc698e694	0.0025	32

## Global dereplication, clustering and chimera detection

Next we pool the 81 samples and do the dereplication, clustering and chimera detection:

```sh
VSEARCH=$(which vsearch)
SWARM=$(which swarm)
TMP_FASTA=$(mktemp --tmpdir=".")
FINAL_FASTA="03_OwnSamples.fas"

# Pool sequences
cat *.fas > "${TMP_FASTA}"

# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
             --threads 46 --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

# Clustering
THREADS=46
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
"${SWARM}" \
    -d 1 -f -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s ${FINAL_FASTA/.fas/_1f.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}

# Sort representatives
"${VSEARCH}" --threads 46 --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/_1f_representatives.fas}
rm ${TMP_REPRESENTATIVES}
  
# Chimera checking
REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
"${VSEARCH}" --threads 46 --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"
```

We get several output files (also see the [Swarm Manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf) for more info):

- `03_OwnSamples.fas`
	- The globally dereplicated .fasta file, containing unique sequences and their counts (again represented by the `size=` in the ID) from all pooled samples
- `03_OwnSamples_1f.stats`
	- Tab-separated statistics file with each OTU in a row with seven columns:
		1. Number of unique amplicons in the OTU
		2. Total abundance of amplicons in the OTU
		3. Seed ID (again that cryptic name)
		4. Initial abundance of that seed
		5. Number of amplicons with an abundance of 1 in the OTU
		6. Maximum number of iterations before the OTU reached its limit
		7. Number of steps along the path joining the seed and the furthermost amplicon in the OTU
- `03_OwnSamples_1f.struct`
	- Internal structure file of all nearly-identical amplicons with five columns:
		1. Amplicon A ID
		2. Amplicon B ID
		3. Number of differences between A and B
		4. OTU number (in order of delineation)
		5. Number of steps from the OTU seed to Amplicon B
- `03_OwnSamples_1f.swarms`
	- List of Amplicons belonging to each OTU (one OTU per line)
- `03_OwnSamples_1f_representatives.fas`
	- Actually the same as `03_OwnSamples.fas`, but with OTUs sorted by abundance
- `03_OwnSamples_1f_representatives.uchime`
	- Chimera detection results with 18 columns (see the [Vsearch Manual](https://manpages.debian.org/stretch/vsearch/vsearch.1.en.html) for details), the last column is the most important for us:
		- Y: query is chimeric
		- N: query is not chimeric
		- ?: Borderline case

## Taxonomic assignment

*To be updated, use this command until a reliable Database is available:*

```sh
grep "^>" 03_OwnSamples_1f_representatives.fas | \
    sed 's/^>//
         s/;size=/\t/
         s/;$/\t100.0\tNA\tNA/' > 03_OwnSamples_1f_representatives.results
```

## Build the OTU Table

The OTU table is built with a python script (`OTU_contingency_table.py`, see below. Note that this is **python 2.7**). All of the previously generated files are used to build the table. Scritly speaking, the resulting table is not a OTU table *per se*, but rather a **contingency table**, because it contains the abundances per OTU and sample (only this would be considered a "real" OTU table) as well as OTU-metadata (like quality, sequence, seed ID and so on).

The files must be provided in a specific order, which is:

- Representatives `03_OwnSamples_1f_representatives.fas`
- Statistics file `03_OwnSamples_1f.stats`
- Swarms `03_OwnSamples_1f.swarms`
- Uchime `03_OwnSamples_1f_representatives.uchime`
- Quality file `01_Oomycota.merged.qual`
- Taxonomic assignment `03_OwnSamples_1f_representatives.results`

So the final command would be:

```sh
python2.7 \
	OTU_contingency_table.py \ 
	03_OwnSamples_1f_representatives.fas \
	03_OwnSamples_1f.stats \ 
	03_OwnSamples_1f.swarms \ 
	03_OwnSamples_1f_representatives.uchime \ 
	01_Oomycota.merged.qual \ 
	03_OwnSamples_1f_representatives.results \
	K*.fas > 04_OwnSamples_OTU_ContingencyTable.tsv
```

`K*.fas` refers to the samples (e.g. K129_Tot.fas etc.), make sure all of them are in the same folder/path.

The `OTU_contingency_table.py` script looks like this:

```python
#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
    Read all fasta files and build a sorted OTU contingency
    table. Usage: python OTU_contingency_table.py [input files]
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr>"
__date__ = "2016/03/07"
__version__ = "$Revision: 5.0"

import os
import re
import sys
import operator

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def representatives_parse():
    """
    Get seed sequences.
    """
    separator = ";size="
    representatives_file = sys.argv[1]
    representatives = dict()
    with open(representatives_file, "rU") as representatives_file:
        for line in representatives_file:
            if line.startswith(">"):
                amplicon = line.strip(">;\n").split(separator)[0]
            else:
                representatives[amplicon] = line.strip()

    return representatives


def stats_parse():
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = sys.argv[2]
    stats = dict()
    seeds = dict()
    with open(stats_file, "rU") as stats_file:
        for line in stats_file:
            cloud, mass, seed, seed_abundance = line.strip().split(separator)[0:4]
            stats[seed] = int(mass)
            seeds[seed] = (int(seed_abundance), int(cloud))
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.iteritems(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats, seeds


def swarms_parse():
    """
    Map OTUs.
    """
    separator = "_[0-9]+|;size=[0-9]+;?| "  # parsing of abundance annotations
    swarms_file = sys.argv[3]
    swarms = dict()
    with open(swarms_file, "rU") as swarms_file:
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            seed = amplicons[0]
            swarms[seed] = [amplicons]

    return swarms


def uchime_parse():
    """
    Map OTU's chimera status.
    """
    separator = " "
    uchime_file = sys.argv[4]
    uchime = dict()
    with open(uchime_file, "rU") as uchime_file:
        for line in uchime_file:
            OTU = line.strip().split("\t")
            try:
                seed = OTU[1].split(";")[0]
            except IndexError:  # deal with partial line (missing seed)
                continue
            try:
                status = OTU[17]
            except IndexError:  # deal with unfinished chimera detection runs
                status = "NA"
            uchime[seed] = status

    return uchime


def quality_parse():
    """
    List good amplicons.
    """
    quality_file = sys.argv[5]
    quality = dict()
    with open(quality_file, "rU") as quality_file:
        for line in quality_file:
            sha1, qual, length = line.strip().split()
            quality[sha1] = float(qual) / int(length)

    return quality


def stampa_parse():
    """
    Map amplicon ids and taxonomic assignments.
    """
    separator = "\t"
    stampa_file = sys.argv[6]
    stampa = dict()
    with open(stampa_file, "rU") as stampa_file:
        for line in stampa_file:
            amplicon, abundance, identity, taxonomy, references = line.strip().split(separator)
            stampa[amplicon] = (identity, taxonomy, references)

    return stampa


def fasta_parse():
    """
    Map amplicon ids, abundances and samples.
    """
    separator = ";size="
    fasta_files = sys.argv[7:]
    samples = dict()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        sample = os.path.basename(fasta_file)
        sample = os.path.splitext(sample)[0]
        samples[sample] = samples.get(sample, 0) + 1
        with open(fasta_file, "rU") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon, abundance = line.strip(">;\n").split(separator)
                    abundance = int(abundance)
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        # deal with duplicated samples
                        amplicons2samples[amplicon][sample] = amplicons2samples[amplicon].get(sample, 0) + abundance
    # deal with duplicated samples
    duplicates = [sample for sample in samples if samples[sample] > 1]
    if duplicates:
        print("Warning: some samples are duplicated", file=sys.stderr)
        print("\n".join(duplicates), file=sys.stderr)
    samples = sorted(samples.keys())

    return amplicons2samples, samples


def print_table(representatives, stats, sorted_stats,
                swarms, uchime, amplicons2samples,
                samples, quality, seeds, stampa):
    """
    Export results.
    """
    # Print table header
    print("OTU", "total", "cloud",
          "amplicon", "length", "abundance",
          "chimera", "spread", "quality",
          "sequence", "identity", "taxonomy", "references",
          "\t".join(samples),
          sep="\t", file=sys.stdout)

    # Print table content
    i = 1
    for seed, abundance in sorted_stats:
        sequence = representatives[seed]
        occurrences = dict([(sample, 0) for sample in samples])
        for amplicons in swarms[seed]:
            for amplicon in amplicons:
                for sample in samples:
                    occurrences[sample] += amplicons2samples[amplicon].get(sample, 0)
        spread = len([occurrences[sample] for sample in samples if occurrences[sample] > 0])
        sequence_abundance, cloud = seeds[seed]

        # Quality
        if seed in quality:
            high_quality = quality[seed]
        else:
            high_quality = "NA"
        
        # Chimera checking (deal with incomplete cases. Is it useful?)
        if seed in uchime:
            chimera_status = uchime[seed]
        else:
            chimera_status = "NA"

        # Chimera checking (deal with incomplete cases. Is it useful?)
        if seed in stampa:
            identity, taxonomy, references = stampa[seed]
        else:
            identity, taxonomy, references = "NA", "NA", "NA"

        # output
        print(i, abundance, cloud,
              seed, len(sequence), sequence_abundance,
              chimera_status, spread, high_quality, sequence,
              identity, taxonomy, references,
              "\t".join([str(occurrences[sample]) for sample in samples]),
              sep="\t", file=sys.stdout)
        i += 1

    return


def main():
    """
    Read all fasta files and build a sorted OTU contingency table.
    """
    # Parse taxonomic results
    representatives = representatives_parse()

    # Parse stats
    stats, sorted_stats, seeds = stats_parse()

    # Parse swarms
    swarms = swarms_parse()

    # Parse uchime
    uchime = uchime_parse()

    # Parse quality
    quality = quality_parse()

    # Parse taxonomic assignment results
    stampa = stampa_parse()
    
    # Parse fasta files
    amplicons2samples, samples = fasta_parse()

    # Print table header
    print_table(representatives, stats, sorted_stats, swarms,
                uchime, amplicons2samples, samples, quality,
                seeds, stampa)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
```

## Filter OTU Table

The OTU table has several columns: 

1.  OTU (each OTU is numbered)
2.  total (total number of reads in the OTU)
3.  cloud (total number of unique sequences in the OTU)
4.  amplicon (identifier of the OTU representative)
5.  length (length of the OTU representative)
6.  abundance (abundance of the OTU representative)
7.  chimera (is it a chimera? Yes, No, ?)
8.  spread (number of samples where the OTU occurs)
9.  quality (minimum expected error observed for the OTU representative, divided by sequence length)
10.  sequence (sequence of the OTU representative)
11.  identity (maximum similarity of the OTU representative with reference sequences)
12.  taxonomy (taxonomic assignment of the OTU representative)
13.  references (reference sequences closest to the OTU representative)

The trailing columns are the abundances of the OTUs in these samples. 

The only filtering we applied so far is removing unmerged reads, very short reads and reads without tags and primers. What we want to do now is to remove chimeras, low quality OTUs and OTUs represented by less than 0.05% of the total reads (this is a bit harsher than removing the classical singletons). Then, we split the Contingency table into the OTU-table and metadata table. To do so, we can use `awk` and `cut`, as described in this script:

```sh
TABLE="04_OwnSamples_OTU_ContingencyTable.tsv"
FILTERED="${TABLE/.tsv/_filtered.tsv}"

head -n 1 "${TABLE}" > "${FILTERED}"
cat "${TABLE}" | awk '$2 >= 141 && $7 == "N" && $9 <= 0.0002 && $5 >= 150' >> "${FILTERED}"

# Filter:
# $2 OTUs with less than 141 Reads (--> 0.05% of total Reads)
# $7 Chimeras
# $9 OTUs with a quality value of less than 0.0002
# $5 Minimum Amplicon length
```

To check how many reads were processed in total, you can run the `rereplicate` function in `vsearch` and count the number of reads:

```sh
vsearch \
	--rereplicate 03_OwnSamples_1f_representatives.fas \
	--threads 46 \
	--output Rereplicated.fas
	
grep -c ">" Rereplicated.fas 
```

	## 2829155
	
Then we multiply this number by 0.00005, which is in this case 142. This means that we remove OTUs with 141 or less reads to get all the OTUs represented by at least 0.05% of all reads.

Next we need to do the taxonomic annotation to check for contaminants and assign species to our OTUs.