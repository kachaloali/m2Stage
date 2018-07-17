export INFILE=tracheophyta.fasta

lipm_fastasplitter.pl --num_per_slice 5000 --outprefix $INFILE.split --in $INFILE

export DB=/db/biomaj/generic/nr/current/blast/

ls $INFILE.split* | xargs -i echo "blastp -query {} -db $DB/nr -num_alignments 20 -max_hsps 20 -num_threads 20  -evalue 1e-5  -show_gis -outfmt 11 -out {}.blastp.nr.asn" > cmd.blastp
