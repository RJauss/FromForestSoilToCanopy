## Why `vsearch` and not `blast`?

Good question! `blast` has a major drawback. It is - as the name suggests - a **local alignment** search tool. This means that a OTU and a reference sequence can have a percent identity of 100%, even if the sequence of the OTU is much smaller than the reference sequence: It does not penalize *terminal gaps*.

Our clustering algorithm, however, does in fact penalize these gaps and counts them as a mismatch. This does in fact make sense, as our target sequence is flanked by conserved regions (i.e. 18S and 5.8S rDNA). Terminal gaps are real sequence information! Maybe this sketch makes it clearer:

```
                28S        ITS1         5.8S
Sequence 1: ...attatc|AATTGATATCTATCTG|ccgcta...
Sequence 2: ...attatc|AATTGATATCTAT|ccgcta...
```
Now when we perform the PCR and remove the conserved ribosomal sequences:
```
Sequence 1: AATTGATATCTATCTG
Sequence 2: AATTGATATCTAT
```
These sequences will **not** be clustered into one OTU, because they have more than 1 mismatch (= terminal gaps). But because they are 100% identical, they will receive **the same** blast hit. The taxonomic diversity of our OTUs might be inflated.

`vsearch`, however, can penalize these terminal gaps and counts them as mismatches. In the above example, the sequences would not be 100% identical, but 81%. The calculation would be as follows:

![equation](https://latex.codecogs.com/gif.latex?%7Bmatching%5C%20columns%20%5Cover%20alignment%5C%20length%7D%20%3D%20%7B13%20%5Cover%2016%7D%20%3D%2081%2C25%5C%25)

This is why we cut the reference sequences with the same primers as we used in our samples. 

## Running vsearch

Here we need the Oomycete sequences which we extracted before as well as the ITS1 database which we built with our primers. Then we can run this command:

```bash
QUERY="04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences.fas"
DATABASE="Oomycota_NCBI_ITS1_PrimerExtracted_cleanHeader.fas"
THREADS=46
OUTPUT="${QUERY/.fas/_vsearch-ITS1-BestHit.tsv}"

 # search for best hits
vsearch \
    --usearch_global ${QUERY} \
    --threads ${THREADS} \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --blast6out ${OUTPUT} \
    --maxaccepts 0 \
    --maxrejects 0 \
    --maxhits 1 \
    --output_no_hits \
    --db ${DATABASE} \
    --id 0.7 \
    --iddef 1 \
    --query_cov 0
```

vsearch reports OTUs with no hit (or less than 70% identity to a reference) with an asterisk. To make it easier for the following steps, we can change the asterisk to a fake taxonomy like this:

```bash
sed -i \
's/\*/XX123456.1 Order_NoHit|Family_NoHit|Genus_NoHit|Genus_NoHit_Species_NoHit/g' \
04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit.tsv 
```

## Refine Annotation

Now we can refine the annotation based on the percent identity. To do so, run this python script:

```python
with open("04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit.tsv", "r") as vsearchOut:
    with open("04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined.tsv", "w") as edited:                       
        for line in vsearchOut:
            line = line.rstrip()
            line = line.split("\t")
            OTU = line[0]
            Hit = line[1]
            percID = line[2]
            Hit = Hit.split(" ")
            AccNo = Hit[0]
            Taxonomy = Hit[1]
            Taxonomy = Taxonomy.split("|")
            Order = Taxonomy[0]
            Family = Taxonomy[1]
            Genus = Taxonomy[2]
            Species = Taxonomy[3]
            if float(percID) < 80:
                if Order != "Order_NoHit":
                    Species = Order + "_Species_NoHit"
                    Genus = Order + "_Genus_NoHit"
                    Family = Order + "_Family_NoHit"
            elif float(percID) < 90:
                Species = Family + "_Species_NoHit"
                Genus = Family + "_Genus_NoHit"
            elif float(percID) < 95:
                Species = Genus + "_Species_NoHit"
            if "sp." in Species:
                Species = Genus + "_Species_NoHit"
            TaxonomyNew = Order + "|" + Family + "|" + Genus + "|" + Species
            edited.write(str(OTU) + "\t" + TaxonomyNew + "\t" + AccNo + "\t" + str(percID) + "\n")
```

For the statistical analyses, it would be good to have the taxonomy tab-separated instead of a pipe. Running this `sed` command turns the refined annotation table into a nice tsv file:

```bash
sed 's/|/\t/g' \
04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined.tsv > \
04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined_noPipe.tsv 
```

## Paste Taxonomy into Contingency Table

Our Contingency Table still contains OTUs which are no oomycetes, and the IDs are just plain numbers. With the following bash script, we filter our precious oomycete OTUs and also rename the IDs based on the taxonomy:

```bash
head -n 1 04_OwnSamples_OTU_ContingencyTable_filtered.tsv > 05_FinalOTU_ContingencyTable.tsv

 #set new line as separator

IFS=$'\n'
set -f

for line in $(cat 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences_vsearch-ITS1-BestHit_AnnotationRefined.tsv); do
	OTU=$(echo "${line}" | cut -f 1)
	Taxonomy=$(echo "${line}" | cut -f 2) 
	grep -w -m 1 "^${OTU}" 04_OwnSamples_OTU_ContingencyTable_filtered.tsv | sed "s/^${OTU}/X${OTU}_${Taxonomy}/" \
	>> 05_FinalOTU_ContingencyTable.tsv; 
done

 #Build the clean OTU table (abundance only) for QIIME's biom format:

cut -f 1,14- 05_FinalOTU_ContingencyTable.tsv > 05_OwnSamples_OTU_Table.tsv

 #Build the Feature Metadata for Qiime:

cut -f -13 05_FinalOTU_ContingencyTable.tsv > 05_OwnSamples_Feature_Metadata.tsv
```