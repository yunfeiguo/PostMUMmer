#!/usr/bin/env python
'''take a query fasta, and optionally a bed
run lastz with hg38 as reference
if a bed is given, then each sequence in query fasta should be mapped a specific reference sequence
'''
import os
import sys
import logging
import tempfile
import time
logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')
#global variables
cleanQ = []
submitInterval = 3 #number of seconds between qstat
checkInterval = 1 #sec
debug = 1
maxConcurrentLaird = 300
maxConcurrentOther = 200
jobReg = {} #record job queue, status, sentinel file path
chr_dir = "/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/hg38_chr/unmasked"
chrFile = "/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/GCA_000001405.15_GRCh38_full_analysis_set.fna"
qsubLaird = "qsub -V -l walltime=23:59:0 -l nodes=1:ppn=1 -A lc_kw -q laird -l mem=2GB -S /bin/bash"
qsubMain = "qsub -V -l walltime=23:59:0 -l nodes=1:ppn=1 -A lc_kw -l mem=2GB -S /bin/bash"
cwd = os.getcwd()

def safeMkdir(dir):
    if not os.path.isdir(dir):
        try:
            os.mkdir(dir)
        except OSError:
            logging.critical("%s exists" % dir)
def readFastaIdx(idx):
    try:
       fh = open(idx,'r')
    except IOError as e:
       logging.critical("I/O error(%d): %s" % e.errno,e.strerror)
    ctg = {}
    for line in fh:
        f = line.split()
        assert len(f) == 5
        (id,length,offset,nCharLine,nByteLine) = f
        ctg[id] = {
                'length':length,
                'offset':offset,
                'nCharLine':nCharLine,
                'nByteLine':nByteLine,
                }
    fh.close()
    return(ctg)
def genTempScript(cmd):
    assert type(cmd) is list
    fd, tmpFile = tempfile.mkstemp(suffix='pairwiseLastz')
    os.write(fd,'#!/bin/bash\n')
    os.write(fd,'set -e\n') #let shell run the script, exit at first error
    os.write(fd,'set -o pipefail\n') #let shell run the script, exit at first error
    for c in cmd:
        os.write(fd,c+'\n')
    os.close(fd)
    cleanQ.append(tmpFile) #clean at the very end
    return(tmpFile)
def clean():
    for f in cleanQ:
        os.remove(f)
def getJobCount(qName):
    '''
    check all jobs not finished
    remove any if finished
    return count for a particular queue
    '''
    #when all jobs are submitted, wait before check
    if len(jobReg) == len([v for k,v in jobReg.items() if v['submitted']]):
        time.sleep(checkInterval)
        #pass
    else:
        time.sleep(submitInterval)
        #pass
    #make sure we mark a job as finished once it is done
    for i in jobReg:
        if os.path.exists(jobReg[i]['doneFile']):
            jobReg[i]['finished'] = True
    total = len([v for k,v in jobReg.items() if v['submitted'] and v['q']==qName])
    finished = len([v for k,v in jobReg.items() if v['finished'] and v['q']==qName])
    return(total-finished)

def submitLastzJob(result_dir, fa, rID, q, qsub):
    '''
    example lastz command
    REF=/home/rcf-02/yunfeigu/proj_dir/pacbio_reference/data/hg38_chr/unmasked/chr1.fa
    QUERY=/auto/rcf-proj/kw/yunfeigu/pacbio_reference/lastz/hx1f2a1_one_ctg_ref_hg38/000000F.fasta
    lastz $REF $QUERY --notransition --step=20 --ambiguous=iupac --format=maf --gappedthresh=1000000 --identity=90 > one_ctg_on_hg38chr1.maf
    
    steps to submit
    	1. mkdir invidual dir
    	2. extract the query sequence
    	3. generate script
    	4. submit"
    '''
    time.sleep(submitInterval)
    dir = os.path.join(result_dir,rID)
    ref = os.path.join(dir,"ref_%s.fa" % rID)
    stde = os.path.join(dir,"stderr")
    stdo = os.path.join(dir,"stdout")
    individualDoneFile = os.path.join(dir,"%s.done" % rID)
    safeMkdir(dir)
    jobReg[rID] = {
            'doneFile':individualDoneFile,
            'q':q,
            'submitted':False,
            'finished':False,
            }
    script = None
    if not os.path.lexists(ref):
        os.symlink(os.path.join(chr_dir,"%s.fa"%rID),ref)
    script = genTempScript([
        "cd $TMPDIR",
		"lastz %s %s --notransition --step=20 --ambiguous=iupac --format=maf --gappedthresh=1000000 --identity=90 > ref_%s_vs_query.maf" % (ref,fa,rID),
		"find . -name '*.maf' -size +4k | xargs -I{} cp {} %s" % (dir),
		"touch %s" % (individualDoneFile),
        ])
    if script is None:
        logging.critical("Nonetype for script")
    if not (os.path.exists(individualDoneFile) or jobReg[rID]['submitted']): 
        os.system("%s -e %s -o %s %s" % (qsub, stde, stdo, script))
        #print("%s -e %s -o %s %s" % (qsub, stde, stdo, script))
        jobReg[rID]['submitted'] = True
#####################################################################
desc = '''
query_mapping.bed has 6 columns, query_ctg_id, start, end, chr, chr_start, chr_end
chr refers to the chromosome the query_ctg_id is mapped to'''
if len(sys.argv) != 3 and len(sys.argv) != 2:
    sys.exit("Usage: %s <query.fa> [query_mapping.bed]\n" % sys.argv[0])
query = sys.argv[1]
bed = None
if len(sys.argv) == 3:
    bed = sys.argv[2]
result_dir = os.path.join(cwd,"lastz_result")
safeMkdir(result_dir)
if bed is not None:
    sys.exit("function for .bed is not yet implemented")
	#when no BED file, align each query contig to every chr, one at a time
allChr = readFastaIdx(chrFile+".fai")
count = 0
for i in allChr:
    count += 1
    if getJobCount('laird') <= maxConcurrentLaird:
        submitLastzJob(result_dir = result_dir, fa = query, rID = i, q = 'laird', qsub = qsubLaird)
    elif getJobCount('other') <= maxConcurrentOther:
        submitLastzJob(result_dir = result_dir, fa = query, rID = i, q = 'other', qsub = qsubMain)
    else:
        count -= 1
clean()
#logging.info("All done")
