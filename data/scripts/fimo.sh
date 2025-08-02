#!/bin/bash

usage()
{
cat << EOF
usage $0 <options>
A TF-Target Validator based on TFBS (transcription factor binding sites) being present on the promoter region sequences of their putative targets
OPTION
 -h	show this Help message
 -n	path to file containing the Co-expression network (delimiter="\t")
 -p	path to file containing the promoter sequences of all genes
 -t	File with LOC ID one per line;
 -v	minimum pvalue to establish a valid match in FIMO (default = 1.0E-4)
 -o	output directory name (default = TF_Targets_match)
 -m directory where i have my loc.meme files saved
DETAILS
This script receives as input: 1)The co-expression network (3 mandatory columns: LOC1 | LOC2 | edge_score/irp_score); 2)The table with the gene promoter sequences (2 mandatory columns: Gene | Sequence); 3)A list with TF IDs (ex:LOC123).
It creates a output directory with severall output files, one of which is a table with 2 columns: TF | Validated Target
The TF putative targets are the genes that the TF is linked to by an edge in the network, and the validation conducted in this script is a check using FIMO (MEME Suite) of a motif match between the TFBS (transcription factor binding site) and the promoter sequences of its putative targets.
Runned this file using: ./fimo.sh -n ../seidr_output/network.txt -p ../cis_elements/promoter_sequences.csv -t ../connecTF/TF_quercus_suber_loc_list.txt

EOF
}

printf "\n\n --------------------  HELLO AND WELCOME TO THE TF-TARGET VALIDATOR  --------------------- \n\n"

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"


#get options
while getopts "h:n:p:t:v:o:m:" OPTION
do
   case $OPTION in
     h) usage; exit 1;;
     n) network=$OPTARG;;
     p) promoters=$OPTARG;;
     t) tf=$OPTARG;;
     v) pvalue=$OPTARG;;
     o) outdir=$OPTARG;;
     m) meme_dir=$OPTARG;;
     ?) usage; exit;;
    esac
done

#Check that all required arguments were passed

if [[ -z $network ]] || [[ -z $promoters ]] || [[ -z $tf ]]
then
	printf "\n-----------------\n  --  ERROR: options -n ; -p and -t are required! -- \n------------------\n\n"
	usage
	exit 1
fi

if [[ -z $pvalue ]];
then
	printf " | Using the FIMO default pvalue threshold of 1.0E-4                        |\n"
	pvalue="1.0E-4"
else
	pvalue=$pvalue
	printf " | Using the user defined pvalue treshold = $pvalue                            |\n"
fi

if [[ -z $outdir ]];
then
	printf " | Using the default output directory = $BASE_DIR/cis_elements/fimo_out             |\n"
	outdir="$BASE_DIR/cis_elements/fimo_out"
else
	outdir=$outdir
	printf " | Outputting the files to the Out-directory = $outdir                         |\n"
fi

if [[ -z $meme_dir ]];
then
	printf " | Using the default meme file directory = $BASE_DIR/cis_elements/TF_motifs             |\n"
	meme_dir="$BASE_DIR/cis_elements/TF_motifs"
else
	meme_dir=$meme_dir
	printf " | Outputting the files to the Out-directory = $meme_dir                         |\n"
fi

#   ---------------------------------       MAIN CODE     --------------------------------   #
mkdir -p $outdir
while read f
do
    FILE=$meme_dir/$f.meme
    if [[ -f "$FILE" ]]; then
		#This script receives all the previous inputs plus the TFBS matrix to prepare a fasta file ready to enter into FIMO
		#this should work without the full path but I am having dificulties resolving it
        $BASE_DIR/scripts/TFBS_promoter_matcher.py --network $network --promoter_seq $promoters --tf $f
		# Extracts target gene names from the TFBS FASTA file, removes '>' and any following text, 
		# and stores the cleaned gene names in TF_${f}_target_list.txt.
		grep '>' "$meme_dir"/fasta/TFBS_TF_${f}_fasta.fa | sed 's/>//g' | sed 's/ .*//g' > $meme_dir/fasta/TF_${f}_target_list.txt

		#This line runs the FIMO tool from MEME Suite and performs the motif match between the TFBS and the promoter region sequences of its putative targets
		
		echo "Running FIMO on TF ${f} for motif matching"
		fimo --oc $outdir/fimo_${f} --verbosity 1 --thresh $pvalue $meme_dir/$f.meme $BASE_DIR/cis_elements/TF_motifs/fasta/TFBS_TF_${f}_fasta.fa
		
		
		echo "Printing all matching targets of the TF ${f}"
		cat $outdir/fimo_${f}/fimo.gff | awk '{print $1,$6}' > $outdir/motif_match_${f}.txt
		#cis_value extracted from FIMO.gff files
		awk -F' ' '!a[$1]++' $outdir/motif_match_${f}.txt > $outdir/best_motif_match_with_score_${f}.txt
		if [[ -f "$outdir/motif_match_${f}.txt" ]]; then
			rm "$outdir/motif_match_${f}.txt"
		fi
		sed 's/ .*//g' $outdir/best_motif_match_with_score_${f}.txt > $outdir/best_motif_match_without_score_${f}.txt
	else
		printf "${f}	NO_TFBS\n" >> $outdir/TFs_with_NO_TFBS.txt
    fi

	#Starts building the Output File Table based on the FIMO results. YES - match | NO - no match | NO_TFBS - no TFBS found
	FILE=$meme_dir/fasta/TF_${f}_target_list.txt
	if [ -f "$FILE" ]; then	
		while read target
		do
			echo "Target: [$target]"
			if grep -Fxq "$target" $outdir/best_motif_match_without_score_${f}.txt
			then
				grep "$target" $outdir/best_motif_match_with_score_${f}.txt > $outdir/fimo_score_prov.txt
				sed 's/.* //g' $outdir/fimo_score_prov.txt > $outdir/fimo_score.txt
				fimo_score=$(<$outdir/fimo_score.txt)
				rm $outdir/fimo_score_prov.txt
				rm $outdir/fimo_score.txt
				printf "${f}	${target}	YES_${fimo_score}\n" >> "$outdir/ALL_TF_Matches_Output.txt"
			else
				printf "${f}	${target}	NO\n" >> "$outdir/ALL_TF_Matches_Output.txt"
			fi
		done < "$meme_dir/fasta/TF_${f}_target_list.txt"

		#rm $outdir/TF_${f}_target_list.txt
		#rm $outdir/best_motif_match_with_score_${f}.txt
		#rm $outdir/best_motif_match_without_score_${f}.txt
		#rm $outdir/$at_tf.meme
		#rm $outdir/TFBS_TF_${f}_fasta.fa
		#rm -r $outdir/fimo_${f}
	fi
done < $tf