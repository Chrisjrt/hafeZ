hafeZ creates a number of output files in different formats.


If a putative active prophage is found hafeZ produces six main ouputs:

- `hafeZ_all_roi_seqs.fasta` = file containing the DNA sequences of all the regions of interest identified
- `hafeZ_summary_all_rois.tsv`  = file containing a summary of info related to each region of interest
- `hafeZ_pyhmmer_hits.tsv` = file containing a list of all region if interest orfs and a description of the pvogs/phrogs they hit
- `hafeZ_prophages_for_xxx.png` = image of zscores per base within contigs where a region of interest was identified with the region highlights (one file per contig containing a region of interest)
- `hafeZ_orfs_aa_XXX.fasta` = fasta file containing amino acid sequence of each orf within the roi
- `hafeZ_orfs_dna_XXX.fasta` = fasta file containing the dna sequence of each orf within the roi

N.B. if the `-Z` flag is used an additional input, `zscores_for_contigXXX.png`, will also be generated which shows the Z-scores of each contig examined (i.e. if input genome contains 100 contigs there will be 100 zscore .png files output. This can be useful if the user wants to manuallly inspect for any potential rois that may be missed under default paramaters. )

If no putative active prophages are found hafeZ will output only an empty `hafeZ_summary_all_rois.tsv` file. 

Additionally, the `hafez__.log` file will contain the logging output, and `logs` directory will contain the logs for all the subprocesses in hafeZ (minimap2, samtools, blast etc).
