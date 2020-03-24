It might make sense to first blast against all available sequences instead of a specific Oomycete ITS1 database. This way you can check for contaminants and remove non-oomycete hits, which would just be a "No Hit" in a smaller database and could be misinterpreted as novel diversity.

Here, we blast against a local copy of the nt database, which we downloaded via: 

```bash
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"

for file in $(ls *gz); do 
     tar -xzvpf ${file}; 
done

rm *gz
```

It is of course not the most reliable database, but by far the largest.

Getting the sequences from the filtered contingency table can be done by running this command:

```bash
cat 04_OwnSamples_OTU_ContingencyTable_filtered.tsv | \
	awk '{print ">"$1"\n"$10}' \
	> 04_OwnSamples_OTU_ContingencyTable_filtered_sequences.fas
```

`$1` prints the OTU ID while `$10` corresponds to the representative sequence of each OTU. Don't forget to manually remove the first "sequence" of the .fas file, this is just the header of the table.

Blasting is quite straightforward, but here we add a new term in the `outfmt` parameter to obtain the `TaxID` from every hit which helps us annotate the OTUs:

```bash
blastn \
	-query 04_OwnSamples_OTU_ContingencyTable_filtered_sequences.fas \
	-db "/media/robin/Data/01_Phd/Databases/NCBI/nt" \
	-out 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted.tsv \
	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" \
	-num_threads 46
```

Then we sort the output and filter for the best hit:

```bash
export LC_ALL=C LC_LANG=C; 
sort \
	-k1,1g \ #sort by OTU number
	-k12,12gr \ #then by bitscore
	-k11,11g \ #then by e-value
	-k3,3gr \ #then by percent identity
	04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted.tsv > \
	04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted.tsv

 #set the OTU ID as a variable and grep the best hit from the sorted blast output
for OTU in $(cut -f1 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted.tsv | sort -u); do 
	grep -w -m 1 "^${OTU}" 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted.tsv; 
done > 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit.tsv
```

## Annotate the OTUs
### Build the Taxonomy Table

Now that we have a blast-output with the taxonomy ID of every hit, we can annotate the OTUs by searching for the ID in NCBI's Taxonomy database. To do so, I wrote a python script using the `Entrez` Package from Biopython:

```python

 # import required packages
import pandas as pd
from Bio import SeqIO, Entrez

 # tell NCBI who you are
Entrez.email = 'your@mail.de'

 # build an empty database, which we later fill with taxonomic information
 # you can choose any taxonomic level you like, here I use:
 # superkingdom, phylum, class, order, family, genus and species
TaxTable = pd.DataFrame(data = "NA", 
                        index = [],
                        columns=["Accession", "PercentIdentity", "superkingdom", 
                                 "phylum", "class", "order", 
                                 "family", "genus", "species"])

 # Define input and output files:
with open("04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit.tsv", "r") as blastout:
    with open("04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable.tsv", "w") as output:
 # loop over lines in the blastoutput and extract the OTU number, the TaxID of the hit and the acession number
        for line in blastout:
            line = line.rstrip()
            line = line.split("\t")
            OTU = line[0]
            TaxID = line[12]
            Acc = line[1]
            PercIdent = line[2]
 # let Entrez search the Taxonomy database with the provided TaxID
            target_handle = Entrez.efetch(db = "Taxonomy", id = str(TaxID), retmode = "xml")
            target_records = Entrez.read(target_handle)
 # extract the taxonomic levels which are stored in a dictionary and put them into the TaxTable
            for DictElement in target_records[0]["LineageEx"]:
                if DictElement["Rank"] in TaxTable.columns:
                    TaxTable.at[OTU,DictElement["Rank"]] = DictElement["ScientificName"]
 # finally write the accession number, percent identity and species annotation into the table
            TaxTable.at[OTU,"Accession"] = Acc
            TaxTable.at[OTU,"PercentIdentity"] = PercIdent
            TaxTable.at[OTU,"species"] = target_records[0]["ScientificName"]
        TaxTable.to_csv(output, sep = "\t", index=True, header = True, na_rep="NA")
```

### Remove Contaminants

There are in fact some sequences not belonging to the Oomycota. First we need to filter the sorted blast output for Oomycetes, then we can use this table to filter the Contingency Table:

```bash
head -n 1 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable.tsv \
> 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes.tsv

grep "Oomycetes" 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable.tsv \
>> 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes.tsv
```

Now we grep the Oomycete OTUs in the contingency table to get the corresponding sequences:

```bash
for a in $(awk '{print $1}' 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes.tsv); do
     grep -m 1 -w -e "^${a}" 04_OwnSamples_OTU_ContingencyTable_filtered.tsv | \
     awk '{print ">"$1"\n"$10}' \
     >> 04_OwnSamples_OTU_ContingencyTable_filtered_sequences_NCBI-nt_blasted_sorted_BestHit-TaxonomyTable_Oomycetes_sequences.fas; 
done
```

With the resulting fasta file, we can now perform a finer taxonomic annotation based on an oomycete specific ITS1 database, as described in the next section.