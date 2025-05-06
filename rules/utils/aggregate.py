""" script for maple pipeline

Uses output data files from rule mutation_analysis for all files being processed, caggregates consensus sequences, genotypes, and UMI id for each barcode in the demux
"""

import pandas as pd
from Bio import SeqIO
import gzip

if __name__=='__main__':
    
    
    # Open the gzipped file and read its contents
    with gzip.open(snakemake.input.consensus_seqs, "rt") as handle:
        # Parse the FASTA file and convert it to a list of records
        records = list(SeqIO.parse(handle, "fasta"))

    # Create a DataFrame from the records
    seqs-df = pd.DataFrame([
        {
            "Seq_ID": record.id,
            "sequence": str(record.seq)
        }
        for record in records
    ])

    seq-ids = pd.read_csv(snakemake.input.seq_ids, index_col=True)
    genotypes = pd.read_csv(snakemake.input.genotypes, index_col=True)
    filename = [snakemake.input.filename]
    genotypes['filename'] = filename
    all_data = [seq-ids['Seq_ID'], genotypes['NT_substitutions', 'AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous', 'NT_insertions', 'NT_deletions'], consensus-seqs]
    seq-ids['filename'] = filename
    agg_df = merge(seq-ids, seqs-df, on = 'Seq_ID', how= 'inner')
	# Merge with genotypes by filename
	agg_df = merge(agg_df, genotypes, on = 'filename', how= 'inner')
	# Optionally, sort the combined DataFrame if required (e.g., by 'seq-id')
	agg_df = combined_df.sort_values(by='filename')
	
	# Write the combined DataFrame to a new CSV file
	agg_df.to_csv(snakemake.output.agg, index=True)  # write row indices

	
    
    
    
    
    
    
    
    
    
    
    # pivot mean enrichment scores from each tag and relabel them so that they can be merged with the genotypes dataframe sequentially
  
    
    for mean_csv in input_list:
        tag = mean_csv.replace('enrichment/','').replace('_enrichment-scores-mean.csv','')
        mean_enrichment = pd.read_csv(mean_csv, index_col=False)
        sample_label, barcode = list(mean_enrichment.columns)[:2]
        mean_enrichment = mean_enrichment.pivot(index=barcode, columns=sample_label, values='mean_enrichment_score')
        if snakemake.params.filter_missing_replicates:
            mean_enrichment = mean_enrichment.dropna(how='any')
        mean_enrichment.columns = [f'mean_enrichment_score_{tag}_' + str(sample) for sample in mean_enrichment.columns]
        mean_enrichment.reset_index(inplace=True)
        # rename barcode column to match genotypes barcode column
        mean_enrichment.rename(columns={barcode: 'barcode(s)'}, inplace=True)

        # merge genotypes and mean_enrichment
        genotypes = pd.merge(genotypes, mean_enrichment, on='barcode(s)', how='left')
    genotypes.to_csv(snakemake.output.genotypes_enrichment, index=False)