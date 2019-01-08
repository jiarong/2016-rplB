#! /bin/bash
set -e
module load screed
module load NumPy
source /mnt/home/guojiaro/Documents/vEnv/qiime_pip/bin/activate

Seqfile=/mnt/lustre_scratch_2012/tg/g/data/glbrcNew/newData/itag/16SandITSrawData/RAW/SSU/allAfterQC/KBS.forQiime.fa
Ref=/mnt/home/guojiaro/Documents/db/qiimeDB/Silva/rep_set/Silva_108_rep_set.fna
Paramsfile=/mnt/home/guojiaro/Documents/db/qiimeDB/params_dir/open_ref_pick_params_SSU.txt
Cpu=4
Mapfile=/mnt/scratch/tg/g/data/glbrcNew/newData/itag/KBS.qiimemap
Cols=Plant,Time,PlantTime


Outdir=$Seqfile.qiimeout

time python ~/Documents/software/gits/qiime/scripts/pick_open_reference_otus.py -i $Seqfile -o $Outdir -r $Ref -f -p $Paramsfile -aO $Cpu

cd $Outdir

filter_taxa_from_otu_table.py -i otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_table_mc2_w_tax_no_pynast_failures_tax_filtered.biom -n f__Mitochondria,c__Chloroplast,p__Viridiplantae
rm -f biom.sum
biom summarize-table -i otu_table_mc2_w_tax_no_pynast_failures_tax_filtered.biom -o biom.sum

Mincount=$(grep '^ Min:' biom.sum |xargs|cut -f2 -d ' ')
Mincount=${Mincount/.*}
echo "Subsample $Mincount seqs"
core_diversity_analyses.py -o core_diversity_analysis.out -i otu_table_mc2_w_tax_no_pynast_failures_tax_filtered.biom -m $Mapfile -t rep_set.tre -p $Paramsfile -c $Cols -e $Mincount --recover_from_failure -aO $Cpu
