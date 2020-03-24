## Download Oomycete ITS Sequences

NCBI has a huge collection of Oomycete ITS sequences. To download all these sequences, search for `"(ITS OR Internal Transcribed Spacer) AND "oomycetes"[porgn:__txid4762])"`, then hit "Send to: File" and download the .fasta file as well as the complete genebank record (we'll need that later).

## Extract the ITS1 Sequence

We are looking for the ITS1 sequence of Oomycetes, but what we downloaded might be the ITS2 or some parts of the 18/28S rDNA. Also, we need the reference sequences to be **exactly the same length** as our OTUs (so we can penalize terminal gaps at our taxonomic annotation). So what we do is an *in silico* PCR to extract exactly the same fragment as we amplified in our samples, using the same primers. Here we can modify the script from the metabarcoding pipeline, using `cutadapt` and `vsearch`:

```bash
INPUT="Oomycota_NCBI.fasta"
PRIMER_F="GCGGAAGGATCATTACCAC"
PRIMER_R="GCTCGCACAHCGATGAAGA"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F}))  # primer match is 100% of primer length
MIN_R=$(( ${#PRIMER_R}))

 # Define binaries, temporary files and output files
CUTADAPT="$(which cutadapt) -j 46 --discard-untrimmed -e 0 --minimum-length ${MIN_LENGTH}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTA=$(mktemp)

 # Reverse complement fastq file
"${VSEARCH}" --quiet --threads 46 \
	--fastx_revcomp "${INPUT}" \
	--fastaout "${INPUT_REVCOMP}"

LOG="${INPUT/.fasta/.log}"
FINAL_FASTA="${INPUT/.fasta/_ITS1_PrimerExtracted.fas}"

 # Trim tags, forward & reverse primers (search normal and antisens)
cat "${INPUT}" "${INPUT_REVCOMP}" | \
	${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
	${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_F}" - 2>> "${LOG}" > "${TMP_FASTA}"


 # Discard sequences containing Ns, convert to fasta
"${VSEARCH}" \
	--quiet --threads 46 \
	--fastx_filter "${TMP_FASTA}" \
	--fastq_maxns 0 \
	--fastaout "${FINAL_FASTA}" 2>> "${LOG}"


 # Clean
rm -f "${INPUT_REVCOMP}" "${TMP_FASTA}"
``` 

The sequence descriptions of NCBI can be really long and confusing, so the next step is to rename the IDs so we have a clear and nice taxonomy per sequence. This is where the genebank file is important:

## Renaming the sequences
### Building a taxonomy database

I wrote a little python script, which opens the genebank file and extracts the taxonomic information into a nice taxonomy table:

```python
from Bio import SeqIO

with open("Oomycota_NCBI.gb", "r") as gb:
    with open("Oomycota_NCBI_Taxonomy.tsv", "w") as out:
	#parse the sequences:
        for record in SeqIO.parse(gb, "genbank"):
	#exclude sequences without a taxonomy (environmental sequences):
            if "environmental samples" not in record.annotations["taxonomy"]:
		#write the sequence ID and get the species name:
                out.write(record.id + "\t")
                species = record.annotations["source"]
                species = species.split(" ")
		#write the taxonomic levels (order, family, ...):
                for a in record.annotations["taxonomy"][3:]:
                    if "unclassified" not in a:
                        out.write(a + "\t")
		#It somehow stops at the genus level, so we manually add the species name we got a few lines earlier:
                out.write(species[1] + "\n")
```

The Output is a big table with the accession numbers and the corresponding taxonomy. The table needs some manual correction, because not every entry has for example a family. 

I changed missing entries to "NA" and deleted the phylum and kingdom level (because all of them are eukaryotes and oomycetes and so on), so that the table has exactly 5 columns:

1. Accession Number
2. Order
3. Family
4. Genus
5. Species

This table can now be used to rename the reference sequences

### Renaming Reference Sequence

To get clean sequence IDs with a straightforward taxonomy, we can run this script:

```python
from Bio import SeqIO

 #open an empty dictionary where we store the taxonomic information per Accession number
TaxonomyDict = {}

 #open the reference taxonomy table
with open("Oomycota_NCBI_Taxonomy_edited.tsv", "r") as database:
 #loop over every entry and extract the taxonomic information
    for line in database:
        line = line.rstrip()
        line = line.split("\t")
        processID = line[0]
        order = line[1]
 #this try-and-except loop skips entries with missing values
        try:
            family = line[2]
            genus = line[3]
            species = line[4]
 #do some manual correction
            if genus == "Paralagenidium":
                family = "incertae_sedis"
            if family == "Lagenaceae":
                order = "incertae_sedis"
            if any(["uncultured" in species, "aff." in species, "cf." in species, "NA" in species]):
                species = "sp."
            if family == "NA":
                family = order + "_FamilyNoHit"
            if genus == "NA":
                genus = family + "_GenusNoHit"
        except IndexError:
            continue
 #paste the taxonomy and store it in the dictionary
        TaxonomyDict[processID] = " " + order + "|" + family + "|" + genus + "|" + genus + "_" + species

 #open the output file and the reference fasta
with open("Oomycota_NCBI_ITS1_PrimerExtracted_cleanHeader.fas", "w") as output:
    with open("Oomycota_NCBI_ITS1_PrimerExtracted.fas", "r") as fas:
        for record in SeqIO.parse(fas, "fasta"):
            ID = record.id
            ID = ID.split(" ")
 #get the accession number and assign the taxonomy based on the dictionary
            AccNo = ID[0]
            if AccNo in TaxonomyDict:
                record.id = AccNo + " " + TaxonomyDict[AccNo]
                record.name = ""
                record.description = ""
                SeqIO.write(record, output, "fasta")
```

The script relies on the presence of every taxonomic level (order, family, genus, species) in the taxonomy file and is sensitive to empty cells. You can either edit the taxonomy table with "NA" instead of missing values, or you can try to ommit these cells with an `try` and `except IndexError` loop. 

The resulting fasta file serves as a nice reference database for our next step: the taxonomic annotation