#get data
r1=EcMG1_ATGTCA_L001_R1_001.fastq.gz
r2=EcMG1_ATGTCA_L001_R2_001.fastq.gz
ref=EcMG.fna
if [ ! -f $ref ]
then
    curl -O https://s3-eu-west-1.amazonaws.com/nxtrim-examples/bacteria/${ref}
fi
if [ ! -f $r1 ]
then
    curl -O https://s3-eu-west-1.amazonaws.com/nxtrim-examples/bacteria/${r1}
fi
if [ ! -f $r2 ]
then
    curl -O https://s3-eu-west-1.amazonaws.com/nxtrim-examples/bacteria/${r2}
fi


##assemble with velvet
time ../nxtrim -1 $r1 -2 $r2  -O EcMG --aggressive
velveth output_dir 61 -short -fastq.gz EcMG.se.fastq.gz -shortPaired2 -fastq.gz EcMG.pe.fastq.gz -shortPaired3 -fastq.gz EcMG.mp.fastq.gz -shortPaired4 -fastq.gz EcMG.unknown.fastq.gz
velvetg output_dir -exp_cov auto -cov_cutoff auto -shortMatePaired4 yes
python n50.py  output_dir/contigs.fa 500
