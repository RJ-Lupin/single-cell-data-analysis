# Add new gene or exogenous sequence to reference genome for cellranger-count
# refer to the official cellranger document page: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr

# 1. save the sequence as fasta file (here, refer to "exogenous.fa")
$ cat exogenous.fa | grep -v "^>" | tr -d "\n" | wc -c # this will return the length of sequence

# 2. make gtf file using the length of sequence 
echo -e 'exogenous_gene_name\tunknown\texon\t1\t$(length_of_sequence)\t.\t+\t.\tgene_id "exogenous_gene_name"; transcript_id "exogenous_gene_name"; gene_name "exogenous_gene_name"; gene_biotype "protein_coding";' > exogenous.gtf

# 3. Add the fasta and gtf to original reference fasta and gtf files
# here, let's assume that this is for mouse reference genome
cp Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39_edited.fa
cat exogenous.fa >> Mus_musculus.GRCm39_edited.fa

cp Mus_musculus.GRCm39.111.gtf Mus_musculus.GRCm39.111_edited.gtf
cat exogenous.gtf >> Mus_musculus.GRCm39.111_edited.gtf

# 4. Run cellranger mkref with edited fasta and gtf files
$ cellranger mkref --genome=GRCm39_edited \
    --fasta=Mus_musculus.GRCm39_edited.fa \
    --genes=Mus_musculus.GRCm39.111_edited.gtf \
    --nthreads=$(nproc)
    
