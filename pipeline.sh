##1. processing preparation

mkdir raw_cutadapt raw_trimmomatic mapped_reads assembled_reads qc_result rseqc_result
path_raw_read=/home/zhouwanting/shujin/90_sample-Geneplus/raw_read
path_raw_cutadapt=/home/suqiang/90-samples/raw_cutadapt
path_raw_trimmomatic=/home/suqiang/90-samples/raw_trimmomatic
path_qc_result=/home/suqiang/90-samples/qc_result
path_rseqc_result=/home/suqiang/90-samples/rseqc_result
bam_dir=/home/suqiang/90-samples/bam_reads


##2. cutadapt

#center1
for i in $(cat $path_raw_read/samplelist.txt)
do
        cutadapt -j 20 -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGAC -A AAGTCGGATCGTAGCCATGTCGTTC -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq.gz -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq.gz $path_raw_read/${i}_raw_1.fq.gz $path_raw_read/${i}_raw_2.fq.gz
done 
echo cut adapter done

#center2
for i in $(cat $path_raw_read/samplelist.txt)
do
        cutadapt -j 20 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq.gz -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq.gz $path_raw_read/${i}_raw_1.fq.gz $path_raw_read/${i}_raw_2.fq.gz
done
echo cut adapter done

#center3
for i in $(cat $path_raw_read/samplelist.txt)
do
        cutadapt -j 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${path_raw_cutadapt}/${i}_cutadapt_R1.fastq.gz -p ${path_raw_cutadapt}/${i}_cutadapt_R2.fastq.gz $path_raw_read/${i}_raw_1.fq.gz $path_raw_read/${i}_raw_2.fq.gz
done
echo cut adapter done

##3.trimmomatic

for i in $(cat $path_raw_read/samplelist.txt);
do
trimmomatic PE -threads 15 $path_raw_cutadapt/${i}_cutadapt_R1.fastq.gz $path_raw_cutadapt/${i}_cutadapt_R2.fastq.gz -baseout $path_raw_trimmomatic/${i}_cutadapt_trim.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20; done


##4. STAR 20221220
for i in $(cat /home/suqiang/90-samples/raw_trimmomatic/samplelist-left.txt)
do
        STAR --runThreadN 10 --genomeDir /home/suqiang/ref/STAR_index/ --readFilesIn $path_raw_trimmomatic/${i}_cutadapt_trim_1P.fastq.gz $path_raw_trimmomatic/${i}_cutadapt_trim_2P.fastq.gz --sjdbOverhang 149 --outFileNamePrefix ${bam_dir}/${i}- --outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --chimOutType Junctions SeparateSAMold --chimSegmentMin 10 --readFilesCommand zcat --quantMode GeneCounts
done
echo STAR done

##5.samtools

for i in $(cat $path_raw_read/samplelist.txt)
do
        samtools sort -@ 10 -o ${bam_dir}/${i}_cutadapt_trim.bam ${sam_dir}/${i}_cutadapt_trim.sam
done
echo samtools done

##6. samtools index *.bam

for i in $(cat $path_raw_read/samplelist.txt)
do
        samtools index -@ 10 ${bam_dir}/${i}_cutadapt_trim.bam
done

##7. output the k-mer-based read count against transcript
for j in $(cat EMP1-isoforms-all-1.txt)
    do
        for i in $(cat samplelist.txt)
            do
            grep ${j} ./${i}/isoforms.fpkm_tracking >> 90-sample_EMP1_${j}.txt

            done
        done



