#!/bin/bash

# Input GFF file
GFF_FILE="reference_genome/Methylocystis_sp_ref.gff"
OUTPUT_FILE="gff_cds_table.tsv"

# Print header
echo -e "Gene\tLocus_tag\tProtein_ID\tProduct\tGO_Function\tGO_Process\tGO_Component\tStart\tEnd\tStrand" > $OUTPUT_FILE

# Extract CDS lines and parse attributes
awk -F'\t' '$3=="CDS" {
    start=$4
    end=$5
    strand=$7
    # Split attributes into array
    n=split($9, a, ";")
    gene=""; locus_tag=""; protein_id=""; product=""; go_function=""; go_process=""; go_component=""
    for(i=1;i<=n;i++){
        split(a[i], b, "=")
        key=b[1]; val=b[2]
        if(key=="gene") gene=val
        if(key=="locus_tag") locus_tag=val
        if(key=="protein_id") protein_id=val
        if(key=="product") product=val
        if(key=="go_function") go_function=val
        if(key=="go_process") go_process=val
        if(key=="go_component") go_component=val
    }
    # Remove GO "|xxxx||IEA" annotations
    gsub(/\|[0-9]+(\|\|[A-Z]+)?/, "", go_function)
    gsub(/\|[0-9]+(\|\|[A-Z]+)?/, "", go_process)
    gsub(/\|[0-9]+(\|\|[A-Z]+)?/, "", go_component)
    # Print tab-delimited line
    print gene"\t"locus_tag"\t"protein_id"\t"product"\t"go_function"\t"go_process"\t"go_component"\t"start"\t"end"\t"strand
}' $GFF_FILE >> $OUTPUT_FILE

echo "CDS table saved as $OUTPUT_FILE"

