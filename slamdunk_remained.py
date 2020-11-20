#!/usr/bin/env python
#aim: submitting jobs in parallel
#author: Zhen Xie
import glob
import subprocess
import re
import os
import sys
import argparse
sys.path.append("/home/xzhen/pipelines/slamseq")
from encode_lib_common import strip_ext_bam

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Wrapper for BSUB job submission.', description='')
    parser.add_argument('--input-dir', default='', type=str, help= 'Path to input files.')
    parser.add_argument('--bedfile', default='/research/rgs01/project_space/yu3grp/Network_JY/yu3grp/LAP_Green/190553SLAMseq/UTR/3pUTR.bed', type=str, help= 'Path to bed files.')
    parser.add_argument('--reference', default='/research/rgs01/project_space/yu3grp/Network_JY/yu3grp/LAP_Green/GRCm38.p6.genome.fa', type=str, help= 'Path to reference.')
    parser.add_argument('--memory', default='32000', type=str, help= 'Memory requested to run the analysis.')
    parser.add_argument('--queue', default='standard', type=str,help='Queue to submit the job in HPCF (use bqueues to choose).')
    parser.add_argument('--out-dir', default='', type=str,help='Output Directory.')
    args = parser.parse_args()
    return args

args = parse_arguments()
bam = glob.glob(args.input_dir+'/*R1_001.bam')

hpcfsubmit = 'bsub ' + '-R ' + '"rusage[mem=' + args.memory + ']" ' + '-q ' + args.queue + ' < '

def create_job_file_pe(bedfile, bamfile, reference, in_dir, out_dir):
    basename = os.path.basename(strip_ext_bam(bam))
    prefix = os.path.join(out_dir,basename)
    job_header = '#!/bin/bash\n'
    job_header += '#BSUB -P SLAMseq\n'
    job_header += '#BSUB -J {}_SLAMseq\n'
    job_header += '#BSUB -n 4\n'
    job_header += '#BSUB -B xzhen@stjude.org\n'
    job_header += '#BSUB -N xzhen@stjude.org\n'
    job_header = job_header.format(basename)
    module1 = 'module load conda3/5.1.0\n'
    module1 += 'source activate slamseq\n'

    job_body1 = 'slamdunk filter -o {} -b {} -t 1 {}'
    job_body1 = job_body1.format(prefix, bedfile, bamfile)

    job_body2 = 'slamdunk snp -o {} -r {} -t 16 {}'
    job_body2 = job_body2.format(prefix, reference, prefix+'/'+basename+ '_filtered.bam')

    job_body3 = 'slamdunk count -o {} -s {} -r {} -b {} -t 1 {}'
    job_body3 = job_body3.format(prefix, prefix+'/*.vcf',reference, bedfile, prefix+'/'+basename+ '_filtered.bam')

    jobfile = prefix + "_remained.sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+module1+'\n'+job_body1+'\n'+job_body2+'\n'+job_body3+'\n')
    return jobfile
    
def submit_job(jobf):
    os.system('{}'.format(hpcfsubmit) + jobf)

for sample in range(0,len(bam)):
    submit_job(create_job_file_pe(args.bedfile, bam[sample], args.reference, args.input_dir, args.out_dir))
    

#run it via: python /home/xzhen/pipelines/slamseq/slamdunk_remained.py --input-dir interleave --out-dir interleave
