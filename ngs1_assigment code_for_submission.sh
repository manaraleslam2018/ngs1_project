
1- extract fastq from sra data 

fastq-dump --outdir . --gzip --split-files SRR8797509.sra



2- Prepare the data

2a) https://bioinf.shenwei.me/seqkit/usage/
[INFO] sample by proportion

source activate ngs1
zcat SRR8797509_1.fastq.gz | seqkit sample -p 0.02 -o sample_r1.fastq.gz
zcat SRR8797509_2.fastq.gz | seqkit sample -p 0.02 -o sample_r2.fastq.gz


2b)- https://bioinf.shenwei.me/seqkit/usage/#split
Split sequences into 50 parts

seqkit split sample_r1.fastq.gz -p 5
seqkit split sample_r2.fastq.gz -p 5

2c)- shufle data and make 100000 shuffled reads 

seqkit shuffle sample_r1.fastq.gz > shuffled_sample_r1.fastq.gz
seqkit shuffle sample_r2.fastq.gz > shuffled_sample_r2.fastq.gz

2D)- split shuffeled sample 1& 2
seqkit split shuffled_sample_r1.fastq.gz -p 5
seqkit split shuffled_sample_r2.fastq.gz -p 5

3- FASTQ Quality Control
For Sample 1 only, use FASTQC to report the difference between S1_1 and S1_2


3a) Install the software

source activate ngs1
conda install -c bioconda fastqc

3B) Run the FASTQC for each read end

mkdir ~/ngs_assigment/FASTQC_sample1_2 && cd ~/ngs_assigment/FASTQC_sample1_2
cp  ~/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_001.fastq.gz .
cp  ~/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_002.fastq.gz .

for f in ~/ngs_assigment/FASTQC_sample1_2/*.fastq.gz;do fastqc -t 1 -f fastq -noextract $f;done



4- Trimming
Mild Trimming for SX_1. {unshuffled}

mkdir ~/ngs_assigment/trimmed && cd ~/ngs_assigment/trimmed 
for r in 1 2 3 4 5;do >f1="/home/manar/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_00${r}.fastq.gz"; >f2="/home/manar/ngs_assigment/sample_r2.fastq.gz.split/sample_r2.part_00${r}.fastq.gz"; >newf1="/home/manar/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_00${r}.pe.trim.fastq.gz"; >newf2="/home/manar/ngs_assigment/sample_r2.fastq.gz.split/sample_r2.part_00${r}.pe.trim.fastq.gz"; >newf1U="/home/manar/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_00${r}.se.trim.fastq.gz"; >newf2U="/home/manar/ngs_assigment/sample_r2.fastq.gz.split/sample_r2.part_00${r}.se.trim.fastq.gz"; >adap="/home/manar/anaconda3/envs/ngs1/share/trimmomatic-0.39-0/adapters"; trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:10 MINLEN:36;done

Agrissive Trimming for SX_2. {shuffled}

mkdir ~/ngs_assigment/shuffled_trimmed && cd ~/ngs_assigment/shuffled_trimmed 
for r in 1 2 3 4 5;do 
f1="/home/manar/ngs_assigment/shuffled_sample_r1.fastq.gz.split/shuffled_sample_r1.part_00${r}.fastq.gz"; 
f2="/home/manar/ngs_assigment/shuffled_sample_r2.fastq.gz.split/shuffled_sample_r2.part_00${r}.fastq.gz"; 
newf1="/home/manar/ngs_assigment/shuffled_sample_r1.fastq.gz.split/shuffled_sample_r1.part_00${r}.pe.trim.fastq.gz"; 
newf2="/home/manar/ngs_assigment/shuffled_sample_r2.fastq.gz.split/shuffled_sample_r2.part_00${r}.pe.trim.fastq.gz";
newf1U="/home/manar/ngs_assigment/shuffled_sample_r1.fastq.gz.split/shuffled_sample_r1.part_00${r}.se.trim.fastq.gz"; 
newf2U="/home/manar/ngs_assigment/shuffled_sample_r2.fastq.gz.split/shuffled_sample_r2.part_00${r}.se.trim.fastq.gz"; 
adap="/home/manar/anaconda3/envs/ngs1/share/trimmomatic-0.39-0/adapters";
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:30 MINLEN:36;
done

5- Alignment by BWA for unshuffled data
5a): Download reference

mkdir ~/ngs_assigment/sample_data && cd ~/ngs_assigment/sample_data
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
gunzip gencode.v29.pc_transcripts.fa.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz


cd ~/workdir/sample_data/
READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS
    do  
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    done


5B) install bwa

source activate ngs1
conda install -c bioconda bwa 

5C) index your genome

mkdir -p ~/ngs_assigment/bwa_align/bwaIndex && cd ~/ngs_assigment/bwa_align/bwaIndex
ln -s ~/ngs_assigment/sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa .
bwa index -a bwtsw gencode.v29.pc_transcripts.chr22.simplified.fa


5D) sequence alignment
cd ~/ngs_assigment/bwa_align/

for r in 1 2 3 4 5

    do
    R1="/home/manar/ngs_assigment/sample_r1.fastq.gz.split/sample_r1.part_00${r}.pe.trim.fastq.gz"
    R2="/home/manar/ngs_assigment/sample_r2.fastq.gz.split/sample_r2.part_00${r}.pe.trim.fastq.gz"
    /usr/bin/time -v bwa mem bwaIndex/gencode.v29.pc_transcripts.chr22.simplified.fa $R1 $R2 >sample_part_00${r}.sam
    done

for r in 1 2 3 4 5
do
    samtools flagstat /home/manar/ngs_assigment/bwa_align/sample_part_00${r}.sam > sample_part_00${r}_stats.out
done

6) HISAT allignment 

6A) installhisat2
conda install -c bioconda hisat2 

6B) index your genome
mkdir -p ~/ngs_assigment/hisat_align/hisatIndex && cd ~/ngs_assigment/hisat_align/hisatIndex
ln -s ~/ngs_assigment/sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa .
ln -s ~/ngs_assigment/sample_data/gencode.v29.annotation.gtf .

hisat2_extract_splice_sites.py ~/ngs_assigment/sample_data/gencode.v29.annotation.gtf > splicesites.tsv
hisat2_extract_exons.py ~/ngs_assigment/sample_data/gencode.v29.annotation.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv gencode.v29.pc_transcripts.chr22.simplified.fa gencode.v29.pc_transcripts.chr22.simplified_indexed


6C) sequence alignment

cd ~/ngs_assigment/hisat_align

for r in 1 2 3 4 5
do
R1="/home/manar/ngs_assigment/shuffled_sample_r1.fastq.gz.split/shuffled_sample_r1.part_00${r}.pe.trim.fastq.gz"
R2="/home/manar/ngs_assigment/shuffled_sample_r2.fastq.gz.split/shuffled_sample_r2.part_00${r}.pe.trim.fastq.gz"
hisat2 -p 1 -x hisatIndex/gencode.v29.pc_transcripts.chr22.simplified_indexed --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuffeled_hisat_aligment${r}.sam
done

for r in 1 2 3 4 5
do
    samtools flagstat /home/manar/ngs_assigment/hisat_align/shuffeled_hisat_aligment${r}.sam > sample_part_00${r}_stats.out
done


6- Assembly
Prepare the SAM file for assembly
conda install samtools

Convert the SAM files from bwa aligment into BAM file 

for r in 1 2 3 4 5
do
samtools view -bS /home/manar/ngs_assigment/bwa_align/sample_part_00${r}.sam > /home/manar/ngs_assigment/bwa_align/sample_part_bwa_aligment00${r}.bam
samtools sort /home/manar/ngs_assigment/bwa_align/sample_part_bwa_aligment00${r}.bam -o /home/manar/ngs_assigment/bwa_align/sample_part_bwa_aligment_sorted${r}.bam
done

Convert the SAM files from Haisat aligment into BAM file 

for r in 1 2 3 4 5
do
samtools view -bS /home/manar/ngs_assigment/hisat_align/shuffeled_hisat_aligment${r}.sam > /home/manar/ngs_assigment/hisat_align/shuffeled_hisat_aligment${r}.bam
samtools sort /home/manar/ngs_assigment/hisat_align/shuffeled_hisat_aligment${r}.bam -o /home/manar/ngs_assigment/hisat_align/shuffeled_hisat_aligment_sorted${r}.bam
done

Assembly for 5 samples shuffled.
