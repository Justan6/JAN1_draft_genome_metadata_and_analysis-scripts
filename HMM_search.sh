### HMMER
# HMMER is often used together with a profile database, such as Pfam or many of the databases that participate in Interpro. But HMMER can also work with query sequences, not just profiles, just like BLAST. For example, you can search a protein query sequence against a database with phmmer, or do an iterative search with jackhmmer. 
# https://mgkit.readthedocs.io/en/0.3.2/pipeline/pipeline-hmmer.html#digital-normalisation-and-qc
# https://jlsteenwyk.com/orthofisher/tutorial/index.html
#http://www.math.chalmers.se/Stat/Bioinfo/Master/Courses/BioinformaticsII/HMM_lab.html


#  Make sure your FASTA file is a multiple sequence alignment (MSA)
mafft --anysymbol NiFe_hydrogenase.fasta > NiFe_hydrogenase.aln.fasta

# Strategy: FFT-NS-2 (Fast but rough) Progressive method (guide trees were built 2 times.)

#Build the HMM profile
hmmbuild NiFe.hmm NiFe_hydrogenase.aln.fasta

#Search for NiFe hits in annotated protein set
hmmsearch -E 1e-50 --tblout NiFe.txt NiFe.hmm ../DRAM_anno/genes.faa
hmmsearch --tblout NiFe.txt NiFe.hmm ../DRAM_anno/genes.faa

#Run hmmpress on your HMM profile and then run hmmscan
hmmpress NiFe.hmm  # Only needed if scanning multiple profiles with hmmscan
hmmscan --tblout scan_results.tbl NiFe.hmm ../DRAM_anno/genes.faa



  
# Retrieve the sequences from the file hmmsearch.output: http://cryptogenomicon.org/extracting-hmmer-results-to-sequence-files-easel-miniapplications.html
# One way to extract the hits is with the -A option to the search programs (hmmsearch, phmmer, jackhmmer). This saves all the domains (subsequences) that were significant (passed the "inclusion thresholds") as one alignment, in Stockholm format. Then we can use esl-reformat to convert the alignment to FASTA format:   

    
esl-sfetch --index  all.genomes.faa 

    hmmsearch -E 1e-50 -A fdhF.sto fdhF.hmm all.genomes.faa 

    esl-reformat fasta fdhF.sto > fdhF.fa





mafft --anysymbol All_NiFe_hydrogenase.fasta > All_NiFe_hydrogenase.aln.fasta

iqtree -s All_NiFe_hydrogenase.aln.fasta -o MAG0 -nt 10 -m MFP -bb 1000



    


