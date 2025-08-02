#!/bin/bash

# Used this comman d to generate the global network: ./seidr.sh -n 35 -p 1.64 -d SLOW
# Used this command to generate the targeted network: ./seidr.sh -n 35 -p 1.64 -d SLOW -t tf_and_tr_file.txt 

# To run seidr.sh, you need to have the following files in the same output directory:
  # 1. Raw_Counts.tsv: This file contains the raw counts of the genes.
  # 2. genes.txt: This file contains the list of genes.
  # 3. tf_and_tr_file if running seidr in targeted mode: This file contains the list of transcription factors and target genes.

BIN_DIR="$(cd "$(dirname "$0")/../../bin" && pwd)"
export PATH="$BIN_DIR/seidr/build:$PATH"
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPTS_DIR="$BASE_DIR/scripts"

usage() {
cat << EOF

Usage: $0 [options]

A Crowd Network generator for visualization of Co-Expression Networks.

Options:
  -h      Show this help message and exit.
  -f      Path to the file containing the Raw_Counts.txt table (delimiter: " ").
  -n      Number of threads to use for the Seidr run.
  -a      Aggregation method for crowd network (default: "irp").
  -o      Output directory for generated files and name of the output network.txt file (default: "seidr_output").
  -p      P-value threshold (default: 1.28, approx. 10%). For 5% use 1.64, and for 1% use 2.32.
  -d      Depth level for network generation. Higher depth increases runtime but extracts more information (default: "VERY_SLOW").
  -t      Path to a file containing a list of target genes (one gene per row). Used for generating a Targeted Network.
  -c      Cutoff value for the crowd network.
  -s      Generate additional statistics files:
            - network_stats.txt (overall network statistics)
            - network_nodestats.txt (node-level statistics)
  -e      Output an additional file with only the edge correlation scores.

Extra Information:
    There are currently four methods of aggregation implemented: (-m borda, -m top, -m top2, -m irp(default))
    There are 4 defined depths available:
      -FAST -> PEARSON SPEARMAN PCOR 
      -MEDIUM -> RAW CLR ARACNE
      -SLOW -> NARROMI PLSNET LLR-ENSEMBLE SVM-ENSEMBLE GENIE3 TIGRESS
      -VERY_SLOW -> EL-ENSEMBLE

Examples:
    ./seidr.sh -n *nºTrheads* [optional options :)]
EOF
}

printf "\n\n ----------------------------- STARTING THE SEIDR SCRIPT ------------------------------- \n\n"

while getopts "h:f:n:a:o:p:d:t:c:s:e:" OPTION
do
  case $OPTION in
    h) usage; exit 1;;
    f) rawcountsfile=$OPTARG;;
    n) thread=$OPTARG;;
    a) aggregate=$OPTARG;;
    o) outdir=$OPTARG;;
    p) pvalue=$OPTARG;;
    d) depth=$OPTARG;;
    t) target=$OPTARG;;
    c) cutoff=$OPTARG;;
    s) stats=$OPTARG;;
    e) edgecorscores=$OPTARG;;
    ?) usage; exit;;
  esac
done

# Function to set default values for optional parameters
set_default_values() {
  # Set default values if variables are not set
  [[ -z $thread ]] && {
    printf "\n---------------------\n -- ERROR: option -n is required!-- \n--------------------\n\n"
    usage
    exit 1
  }

  rawcountsfile="${rawcountsfile:-$BASE_DIR/Raw_Counts.txt}"
  printf " | Raw_Counts.txt file found at  = $rawcountsfile\n"
  
  # Default aggregation method
  aggregate="${aggregate:-irp}"
  printf " | Using the aggregation method: $aggregate\n"

  # Default output directory
  outdir="${outdir:-$BASE_DIR/seidr_output}"
  mkdir -p $outdir
  printf " | Outputting the files to the Out-Directory = $outdir\n"

  # Default p-value
  pvalue="${pvalue:-1.28}"
  printf " | Applying p-value for pruning noise: $pvalue\n"

  # Default depth setting
  depth="${depth:-VERY_SLOW}"
  printf " | Applying computational speed: $depth\n"

  # Default target setting
  target="${target:-FALSE}"
  printf " | Targeted network: $target\n"

  # Default cutoff setting
  cutoff="${cutoff:-None}"
  printf " | Applying cutoff: $cutoff\n"

  # Default stats setting
  stats="${stats:-TRUE}"
  if [[ "$stats" != "TRUE" && "$stats" != "FALSE" ]]; then
    printf " | ERROR: Invalid value for -s/--stats. Acceptable values are 'TRUE' or 'FALSE'. Using default: TRUE.\n"
    stats="TRUE"  # Default to TRUE if invalid value
  fi
  printf " | Calculate network and node statistics: $stats\n"

  # Default edge correlation scores
  edgecorscores="${edgecorscores:-TRUE}"
  if [[ "$edgecorscores" != "TRUE" && "$edgecorscores" != "FALSE" ]]; then
    printf " | ERROR: Invalid value for -e/--edgecorscores. Acceptable values are 'TRUE' or 'FALSE'. Using default: TRUE.\n"
    edgecorscores="TRUE"  # Default to TRUE if invalid value
  fi
  printf " | Generate edge correlation scores: $edgecorscores\n"
}

set_default_values

# --------------------------------------------------------- MAIN CODE -------------------------------------------------------

mkdir -p $outdir

#This script receives a Raw_Counts.txt and outputs an expression.tsv and a genes.txt file to the specified outdir.
#These two files are the base material for all calculus for the Seidr crowd network generation toolkit.

#raw_counts_normalization.R --rawcountsfile $rawcountsfile --output $outdir


#for more information -> https://seidr.readthedocs.io/en/latest/source/getting_started/getting_started.html

printf "\n\n | --------------------- Started computation with $thread threads and a $depth depth ---------------------- |\n\n"

#Checks if the network to be generated is in targeted mode or not
if [[ $target == 'FALSE' ]]
then


# FAST
printf "Calculating the Pearson Correlation scores 'fast'.\n"
correlation -m pearson -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/pearson_scores.tsv
seidr import -A -r -u -n PEARSON -o $outdir/pearson_scores.sf -F lm -i $outdir/pearson_scores.tsv -g $outdir/genes.txt
printf "Calculating the Spearman Correlation scores 'fast'.\n"
correlation -m spearman -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/spearman_scores.tsv
seidr import -A -r -u -n SPEARMAN -o $outdir/spearman_scores.sf -F lm -i $outdir/spearman_scores.tsv -g $outdir/genes.txt
printf "Calculating the PCorrelation scores 'fast'.\n"
pcor -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/pcor_scores.tsv
seidr import -A -r -u -n PCOR -o $outdir/pcor_scores.sf -F lm -i $outdir/pcor_scores.tsv -g $outdir/genes.txt

#MEDIUM
if [[ $depth == 'MEDIUM' ]] || [[ $depth == 'SLOW' ]] || [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the RAW scores 'medium'.\n"
mi -m RAW -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/mi_scores.tsv
seidr import -r -u -n MI -o $outdir/mi_scores.sf -F lm -i $outdir/mi_scores.tsv -g $outdir/genes.txt
printf "Calculating the CLR scores 'medium'.\n"
mi -m CLR -i $outdir/expression.tsv -g $outdir/genes.txt -M $outdir/mi_scores.tsv -o $outdir/clr_scores.tsv
seidr import -r -u -z -n CLR -o $outdir/clr_scores.sf -F lm -i $outdir/clr_scores.tsv -g $outdir/genes.txt
printf "Calculating the ARACNE scores 'medium'.\n"
mi -m ARACNE -i $outdir/expression.tsv -g $outdir/genes.txt -M $outdir/mi_scores.tsv -o $outdir/aracne_scores.tsv
seidr import -r -u -z -n ARACNE -o $outdir/aracne_scores.sf -F lm -i $outdir/aracne_scores.tsv -g $outdir/genes.txt
fi

#SLOW
if [[ $depth == 'SLOW' ]] || [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the Narromi scores 'slow'\n"
narromi -O $thread -m interior-point -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/narromi_scores.tsv
seidr import -r -z -n NARROMI -o $outdir/narromi_scores.sf -F m -i $outdir/narromi_scores.tsv -g $outdir/genes.txt
printf "Calculating the Plsnet scores 'slow'\n"
plsnet -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/plsnet_scores.tsv
seidr import -r -z -n PLSNET -o $outdir/plsnet_scores.sf -F m -i $outdir/plsnet_scores.tsv -g $outdir/genes.txt
printf "Calculating the LLR-Ensemble scores 'slow'\n"
llr-ensemble -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/llr_scores.tsv
seidr import -r -z -n LLR -o $outdir/llr_scores.sf -F m -i $outdir/llr_scores.tsv -g $outdir/genes.txt
printf "Calculating the SVM-Ensemble scores 'slow'\n"
svm-ensemble -O $thread -k POLY -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/svm_scores.tsv
seidr import -r -z -n SVM -o $outdir/svm_scores.sf -F m -i $outdir/svm_scores.tsv -g $outdir/genes.txt
printf "Calculating the Genie3 scores 'slow'.\n"
genie3 -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/genie3_scores.tsv
seidr import -r -z -n GENIE3 -o $outdir/genie3_scores.sf -F m -i $outdir/genie3_scores.tsv -g $outdir/genes.txt
printf "Calculating the Tigress scores 'slow'\n"
tigress -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/tigress_scores.tsv
seidr import -r -z -n TIGRESS -o $outdir/tigress_scores.sf -F m -i $outdir/tigress_scores.tsv -g $outdir/genes.txt
fi

#VERY_SLOW
if [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the El-Ensemble scores 'very slow'.\n"
el-ensemble -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/elnet_scores.tsv
seidr import -r -z -n ELNET -o $outdir/elnet_scores.sf -F m -i $outdir/elnet_scores.tsv -g $outdir/genes.txt
fi

else

#-------------------------------------- Running the Algorithms in targeted mode --------------------------------------#

# FAST
printf "Calculating the Pearson Correlation scores 'fast'.\n"
correlation -t $target -m pearson -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/pearson_scores.tsv
seidr import -A -r -u -n PEARSON -o $outdir/pearson_scores.sf -F el -i $outdir/pearson_scores.tsv -g $outdir/genes.txt
printf "Calculating the Spearman Correlation scores 'fast'.\n"
correlation -t $target -m spearman -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/spearman_scores.tsv
seidr import -A -r -u -n SPEARMAN -o $outdir/spearman_scores.sf -F el -i $outdir/spearman_scores.tsv -g $outdir/genes.txt
printf "Calculating the PCorrelation scores 'fast'.\n"
pcor -t $target -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/pcor_scores.tsv
seidr import -A -r -u -n PCOR -o $outdir/pcor_scores.sf -F el -i $outdir/pcor_scores.tsv -g $outdir/genes.txt

#MEDIUM
if [[ $depth == 'MEDIUM' ]] || [[ $depth == 'SLOW' ]] || [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the RAW scores 'medium'.\n"
mi -t $target -m RAW -i $outdir/expression.tsv -g $outdir/genes.txt -M $outdir/mi_full_scores.tsv -o $outdir/mi_scores.tsv
seidr import -r -u -n MI -o $outdir/mi_scores.sf -F el -i $outdir/mi_scores.tsv -g $outdir/genes.txt
printf "Calculating the CLR scores 'medium'.\n"
mi -t $target -m CLR -i $outdir/expression.tsv -g $outdir/genes.txt -M $outdir/mi_full_scores.tsv -o $outdir/clr_scores.tsv
seidr import -r -u -z -n CLR -o $outdir/clr_scores.sf -F el -i $outdir/clr_scores.tsv -g $outdir/genes.txt
printf "Calculating the ARACNE scores 'medium'.\n"
mi -t $target -m ARACNE -i $outdir/expression.tsv -g $outdir/genes.txt -M $outdir/mi_full_scores.tsv -o $outdir/aracne_scores.tsv
seidr import -r -u -z -n ARACNE -o $outdir/aracne_scores.sf -F el -i $outdir/aracne_scores.tsv -g $outdir/genes.txt
fi

#SLOW
if [[ $depth == 'SLOW' ]] || [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the Narromi scores 'slow'\n"
narromi -t $target -O $thread -m interior-point -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/narromi_scores.tsv
seidr import -r -z -n NARROMI -o $outdir/narromi_scores.sf -F el -i $outdir/narromi_scores.tsv -g $outdir/genes.txt
printf "Calculating the Plsnet scores 'slow'\n"
plsnet -t $target -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/plsnet_scores.tsv
seidr import -r -z -n PLSNET -o $outdir/plsnet_scores.sf -F el -i $outdir/plsnet_scores.tsv -g $outdir/genes.txt
printf "Calculating the LLR-Ensemble scores 'slow'\n"
llr-ensemble -t $target -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/llr_scores.tsv
seidr import -r -z -n LLR -o $outdir/llr_scores.sf -F el -i $outdir/llr_scores.tsv -g $outdir/genes.txt
printf "Calculating the SVM-Ensemble scores 'slow'\n"
svm-ensemble -t $target -O $thread -k POLY -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/svm_scores.tsv
seidr import -r -z -n SVM -o $outdir/svm_scores.sf -F el -i $outdir/svm_scores.tsv -g $outdir/genes.txt
printf "Calculating the Genie3 scores 'slow'.\n"
genie3 -t $target -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/genie3_scores.tsv
seidr import -r -z -n GENIE3 -o $outdir/genie3_scores.sf -F el -i $outdir/genie3_scores.tsv -g $outdir/genes.txt
printf "Calculating the Tigress scores 'slow'\n"
tigress -t $target -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/tigress_scores.tsv
seidr import -r -z -n TIGRESS -o $outdir/tigress_scores.sf -F el -i $outdir/tigress_scores.tsv -g $outdir/genes.txt
fi

#VERY_SLOW
if [[ $depth == 'VERY_SLOW' ]]
then

printf "Calculating the El-Ensemble scores 'very slow'.\n"
el-ensemble -t $target -O $thread -i $outdir/expression.tsv -g $outdir/genes.txt -o $outdir/elnet_scores.tsv
seidr import -r -z -n ELNET -o $outdir/elnet_scores.sf -F el -i $outdir/elnet_scores.tsv -g $outdir/genes.txt
fi

fi

#Aggregating the networks

#FAST
echo "Agregating algorithms"

if [[ $depth == 'FAST' ]]; then seidr aggregate -o $outdir/network.sf -m $aggregate $outdir/pcor_scores.sf $outdir/pearson_scores.sf $outdir/spearman_scores.sf; fi

#MEDIUM
if [[ $depth == 'MEDIUM' ]]; then seidr aggregate -o $outdir/network.sf -m $aggregate $outdir/aracne_scores.sf $outdir/clr_scores.sf $outdir/mi_scores.sf $outdir/pcor_scores.sf $outdir/pearson_scores.sf $outdir/spearman_scores.sf; fi

#SLOW
if [[ $depth == 'SLOW' ]]; then seidr aggregate -o $outdir/network.sf -m $aggregate $outdir/aracne_scores.sf $outdir/clr_scores.sf $outdir/genie3_scores.sf $outdir/llr_scores.sf $outdir/mi_scores.sf $outdir/narromi_scores.sf $outdir/pcor_scores.sf $outdir/pearson_scores.sf $outdir/plsnet_scores.sf $outdir/spearman_scores.sf $outdir/svm_scores.sf $outdir/tigress_scores.sf; fi

#VERY SLOW
if [[ $depth == 'VERY_SLOW' ]]; then seidr aggregate -o $outdir/network.sf -m $aggregate $outdir/aracne_scores.sf $outdir/clr_scores.sf $outdir/elnet_scores.sf $outdir/genie3_scores.sf $outdir/llr_scores.sf $outdir/mi_scores.sf $outdir/narromi_scores.sf $outdir/pcor_scores.sf $outdir/pearson_scores.sf $outdir/plsnet_scores.sf $outdir/spearman_scores.sf $outdir/svm_scores.sf $outdir/tigress_scores.sf; fi

#----------------------------------------------------- POSTPROCESSING ----------------------------------------------------

#Pruning noise from the network, dropping edges bellow a specific P-value (1.28 or 0.10% by default)
#If the desired P-value is 0.05% the -p should be 1.64, and for 0.01% the -p should be 2.32
seidr backbone -F $pvalue $outdir/network.sf

seidr view --column-headers -o $outdir/cyto1_network.txt $outdir/network.bb.sf
#replace all semicolons (;) with tab characters in the file
sed 's/;/\t/g' $outdir/cyto1_network.txt > $outdir/cyto2_network.txt

#To discard the NC_Scores use:
#awk '{print $1,$2,$(NF-3)}' cyto2_network.txt > cyto3_network.txt
#To keep the NC_Scores use:
awk '{print $1,$2,$(NF-3),$(NF-1),$(NF)}' $outdir/cyto2_network.txt > $outdir/network.txt

#Applying a cutoff to the network score with network_cut.py if -cutoff value is specified
if [[ $cutoff != 'None' ]]; then python3 $SCRIPTS_DIR/network_cut.py --network $outdir/network.txt --cutoff $cutoff --method ${aggregate}_score --out $outdir/network_${cutoff}.txt; fi

#Option to keep the graph and node statistics into a *$outdir*_network_stats.txt and *$outdir*_network_nodestats.txt
if [[ $stats == 'TRUE' ]]
then
	
#seidr reheader $PWD/$outdir/network.bb.sf
#seidr graphstats -o $PWD/$outdir/${outdir}_network_stats.txt $PWD/$outdir/network.bb.sf
seidr stats --metrics PR,BTW,CLO $outdir/network.bb.sf
seidr view --centrality -c -o $outdir/network_nodestats.txt $outdir/network.bb.sf
fi

#Option to collect the correlation scores in a separate file called *name*_network_corscores.txt.
#This option can be used to import the scores in Cytoscape and color the edges of the network, giving a sense of positive and negative correlation of the network.
if [[ $edgecorscores == 'TRUE' ]] 
then

if [[ $depth == 'FAST' ]]
then
        awk '{print $1,$2,$4,$6,$8}' $outdir/cyto2_network.txt > $outdir/cor1_scores.txt
	#awk '{print $1,$2,$4,$5,$6 }' cyto2_network.txt > cor1_scores.txt
fi

if [[ $depth == 'MEDIUM' ]]
then
	awk '{print $1,$2,$7,$8,$9}' $outdir/cyto2_network.txt > $outdir/cor1_scores.txt
fi

if [[ $depth == "SLOW" ]]
then
	awk '{print $1,$2,$16,$18,$22}' $outdir/cyto2_network.txt > $outdir/cor1_scores.txt
	#awk '{print $1,$2,$10,$11,$14 }' cyto2_network.txt > cor1_scores.txt
fi

if [[ $depth == "VERY_SLOW" ]]
then
	awk '{print $1,$2,$11,$12,$14 }' $outdir/cyto2_network.txt > $outdir/cor1_scores.txt
fi

sed 's/  */\t/g' $outdir/cor1_scores.txt > $outdir/cor2_scores.txt

#Only needed to make sure that the correlations are all positive or all negative for the means
sed "s/;[[:digit:]][[:digit:]]*\.[[:digit:]]//g" $outdir/cor1_scores.txt > $outdir/cor2_scores.txt
sed "s/;[[:digit:]][[:digit:]]*//g" $outdir/cor2_scores.txt > $outdir/cor3_scores.txt
echo "Source Target PCOR_score PEARSON_score SPEARMAN_score" > $outdir/cor4_scores.txt
grep -e "[[:space:]]0.*[[:space:]]0.*[[:space:]]0" $outdir/cor3_scores.txt >> $outdir/cor4_scores.txt
grep -e "-0.*-0.*-0" $outdir/cor3_scores.txt >> $outdir/cor4_scores.txt

#This script calculates the means of the correlation scores in cor1_scores and makes an extra key collumn to allow the import in Cytoscape
python3 $SCRIPTS_DIR/network_scores.py --cor_scores $outdir/cor2_scores.txt --out $outdir/network_corscores.txt
fi

sed 's/  */\t/g' $outdir/network.txt > $outdir/network.csv

#Lastly, eliminate intermediate files.
#rm cyto1_network.txt
#rm cyto2_network.txt
#rm cor1_scores.txt
#rm cor2_scores.txt


#This python script applies a cutoff to the network by their correlation scores if specified
if [[ "$cutoff"  != "None" ]]; then python3 $SCRIPTS_DIR/network_cut.py --network network_corscores.txt --cutoff $cutoff --out cork_cor_scores_mean_cut.txt; else echo "Applying no cutoff to the correlation scores"; fi

printf "\n\n | -------------------------- END OF SEIDR_SCRIPT ----------------------------- |\n"
#-------------------------------------------------- END OF MAIN CODE -----------------------------------------------------

