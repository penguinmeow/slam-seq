# slam-seq
Using SLAMDUNK analyze paired-end slam-set data

#ngm_bsub.py
manage to submit jobs in paralell to mapp the forward and reverse reads of a paired-end fastq file to the reference genome.

#slamdunk_remained.py
using same scripts to submit jobs but weird Error occurs as followed:

python /home/xzhen/pipelines/slamseq/slamdunk_remained.py --input-dir interleave --out-dir interleave
print(strip_ext_bam(bam))     #['interleave/WTBMDM1_suffled_R1_001.bam']
print(strip_ext_fastq(fq))    #['interleave/1953895_WT_BMDM_1_S110_L004_R1_001.fastq.gz']
print(os.path.basename(strip_ext_bam(bam)))                 #WTBMDM1_suffled_R1_001.bam']
print(out_dir)                #interleave
print(prefix)                 #interleave/WTBMDM1_suffled_R1_001.bam']
print(prefix+'remained.sh')   #interleave/WTBMDM1_suffled_R1_001.bam']remained.sh

#slamdunk_bsub.py
using sam scripts and manage to submit the test file, not sure if it can run smoothly as the jobs is pending.
