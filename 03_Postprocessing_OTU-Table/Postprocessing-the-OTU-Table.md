## Extract Sequences from Filtered Table

After we subsampled our OTU table, some OTUs might have been completely removed. It might be useful to get an updated .fasta file of all remaining representative sequences.

Here I wrote a workaround to extract sequences from the subsampled OTU table. This is a bit tricky, because we have to find the remaining OTUs in the contingency- or feature-metadata table. There really has to be a simpler way than mine, but nevertheless it works for me:

```bash
for a in $(awk '{print $1}' 05_OwnSamples_OTU_Table_min-freq-9588.tsv); do
     grep -m 1 -w -e "^${a}" 04_OwnSamples_Feature_Metadata.tsv | \
     awk '{print ">"$1"\n"$10}' \
     >> 05_OwnSamples_OTU_Table_min-freq-9588.sequences.fas; 
done
```

What this script does is to print the OTU number of the filtered table (the number remains the same after filtering, which is good), then we `grep` this number in the beginning of each line in the feature-metadata table (alternatively, the contingency table), pipe it into `awk` and print the ">", the OTU number, a new line and then the sequence - which is in column no. 10 of the metadata table -  and *voilà*, we have a subsampled fasta file.

Maybe not an elegant way, but it is efficient.

## Paste Filtered OTU Table and Filtered Metadata

We subsampled the OTU table to exclude samples with low coverage. Now we also need to subsample the metadata file in order to paste it into the filtered OTU table. Here is how I do it:

```bash
for a in $(csvtool transpose -t TAB -u TAB 05_OwnSamples_OTU_Table_min-freq-9588.tsv | cut -f 1); do
	   grep "${a}" 04_OwnSamples_Sample-Metadata.tsv >> \
	   04_OwnSamples_Sample-Metadata_min-freq-9588.tsv; 
done

csvtool transpose -t TAB -u TAB 05_OwnSamples_OTU_Table_min-freq-9588.tsv | \
cut -f 2- | \
paste 04_OwnSamples_Sample-Metadata_min-freq-9588.tsv - \
> 05_OwnSamples_OTU_Table_min-freq-9588_transposed_withMetadata.tsv
```