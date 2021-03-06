# BrainSpan-analysis

This folder includes code and data for the BrainSpan analysis section of *Effects of 3D culturing conditions on the transcriptomic profile of stem-cell-derived neurons*

# Comparison to BrainSpan dataset 
In order to compare a new RNA-seq to an existing transcriptome profile, such as BrainSpan, ideally the datasets should be processed uniformly. Unfortunately to our knowledge raw data for the BrainSpan transcriptome profiling dataset, in the form of fastq or BAM files are not available for download. Therefore, in an effort to minimize technical differences between the experimental and reference datasets, RNA-seq data of 3D cultures of iN cells was reprocessed to more closely match that of BrainSpan following the alignment and gene quantification protocol described. Ths included data alignment and expression quantifications. Our pipeline was designed to mimic BrainSpan, using their documentation available here: http://help.brain-map.org/download/attachments/3506181/Transcriptome_Profiling.pdf?version=1&modificationDate=1382036562736

## 1) Alignment
We matched BrainSpan by using the same Tophat and gencode versions. Below is a sample alignment script using
-samtools v.0.1.9
-Bowtie v.0.12.9
-Tophat v.2.0.14

`/seq/software/tophat2/tophat  --bowtie1 -r 150 -p 8
-G gencode.v10.fixed.annotation.gtf
./Homo_sapiens_assembly19.fasta ${sample_name}_R1_001.fastq,${sample_name}_R2_001.fastq`

the gencode v10 is available in /data/

## 2) Quantification
Expression in BrainSpan was quantified via RSEQtools, we again followed the pipeline laid out in the BrainSpan documentation above to get as close to their pipeline as possible. This requires converting BAMs into SAMs and then MRFs, running the RSEQtools quantifier and merging the per-sample files. 

a) Turn BAMs into SAMs with SortSam

`samtools sort -o ${sample_name}.sorted.sam $bam`

b) Turn SAM files into MRF files 

`RSEQtools/bin/sam2mrf < $sam > /${sample_name}.sorted.mrf`

c) Run transcript quantification

`RSEQtools/bin/mrfQuantifier gencode.v10.merged.transcripts.txt  multipleOverlap < $mrf > ${sample_name}.expression.txt`

You can then follow sample code in scripts/gather_per_sample_quantifications.R to merge quantifications. The resulting quantifications used in the study are available in data/.
