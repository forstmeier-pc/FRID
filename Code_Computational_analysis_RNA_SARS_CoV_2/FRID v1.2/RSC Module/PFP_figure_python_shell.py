import os
import argparse


#shell to run Rscript to draw PFP figures

def get_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-c', type=str, help='input candidate regions in standard format')
    parser.add_argument('-r', type=str, help='input path to ROI file generated by RSC_part_one.py')
    #parser.add_argument('-ct', type=str, help='input path to ct file of the WT sequence generated by RSC_part_one.py')
    
    args = parser.parse_args()
    
    return args

class Candidate:
    def __init__(self, start, stop, sequence, label):
        self.start = int(start)
        self.stop = int(stop)
        self.label = label
        self.sequence = sequence
        self.length = len(sequence)
        self.list_snps = []

def extract_candidates(PK_candidate_file):
    region_file = open(PK_candidate_file,'r')
    lines = region_file.readlines()
    list_candidates = []
    for line in lines[1:]:
        info = line.split('\t')
        region = Candidate(info[0], info[1], info[2], info[3])
        list_candidates.append(region)
    return list_candidates

def master():
    args = get_args()
    list_cands = extract_candidates(str(args.c))
    for cand in list_cands:
        ct_filename = 'wt_region_'+str(cand.region)+'_'+str(cand.region_stop)+'.ct'
        os.system('Rscript %s %s %s %s' % (int(args.st), int(args.sp), str(args.r),ct_filename))