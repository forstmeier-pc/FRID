#SNPfold_ACI.py
#this script is for ACI usage and serves as a work-around so that the .pbs file can switch python versions while running
#many of functions here are derived from blast_output_reader.ACI
import os
import multiprocessing as mp
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-i', type=str, help='input file')
parser.add_argument('-p', type=str, help='path to SNPfold directory')
parser.add_argument('-c', type=str, help='input candidate regions in standard format')
parser.add_argument('-t',type=int, default=1, help='input the threshold for fraction of pseudoknot base pairs disrupted to be considered a pseudoSNitch. Default = 1.0 (any disruption is a pseudoSNitch)')
parser.add_argument('-st',type=float, default=0.8, help='input the threshold for the correlation coefficient of SNPfold to be considered a riboSNitch. Default = 0.8 (any disruption is a pseudoSNitch)')
parser.add_argument('-rmin',type=int, default=0.1, help='input the threshold for the p-value of the correlation coefficient of RNAsnp to be considered a riboSNitch. Default = 0.1 (any disruption is a pseudoSNitch)')
parser.add_argument('-dmax',type=int, default=0.11, help='input the threshold for p-value of the distance calculated by RNAsnp to be considered a riboSNitch. Default = 0.1 (any disruption is a pseudoSNitch)')

args = parser.parse_args()

input = (args.i)
global path_to_snpfold
path_to_snpfold = args.p

def clean_fasta(SNP):
    output_file_name_clean = 'wt_'+SNP.ID+'_clean.fasta'
    output_file = open(output_file_name_clean, 'w')
    output_file.write(SNP.wt_window)
    output_file.close()
    return output_file_name_clean

def make_snp_file(SNP):
    snp_filename = SNP.ID +'.snp'
    snp_file = open(snp_filename, 'w')
    snp_file.write(SNP.wt+str(61)+SNP.alt)
    snp_file.close()
    return snp_filename

def run_SNPfold(SNP):
    orig_dir = os.getcwd()
    os.chdir(path_to_snpfold)
    snpfold_input_file_name = clean_fasta(SNP)
    SNPfold_output_name = 'SNPfold_output_%s.txt' % (SNP.ID)
    SNPfold_output = open(SNPfold_output_name, 'w')
    snp_file = make_snp_file(SNP)
    python_2_command = 'python2.7 SNPfold_commandline.py %s %s > %s' % (snpfold_input_file_name, snp_file, SNPfold_output_name)
    os.system(python_2_command)
    SNPfold_output = open(SNPfold_output_name, 'r')
    lines = SNPfold_output.readlines()
    if lines[0][0] == 'E':
        SNP.cc_bpprob = ''
        SNP.cc_shannon = ''
        os.remove(SNPfold_output_name)
        os.remove(snpfold_input_file_name)
        os.remove(snp_file)
        os.chdir(orig_dir)
        return SNP
    data = lines[1].split('\t')
    SNP.cc_bpprob = data[1]
    SNP.cc_shannon = data[2].replace('\n','')
    os.remove(SNPfold_output_name)
    os.remove(snpfold_input_file_name)
    os.remove(snp_file)
    os.chdir(orig_dir)
    return SNP

def construct_snp_list(snp_list):
    new_list = []
    temp_item = ''
    for i in snp_list:
        if (i == '[') or (i == ' ') or (i == "'"):
            continue
        elif i == ',':
            new_list.append(temp_item)
            temp_item = ''
        elif i == ']':
            break
        else:
            temp_item += i
    return new_list

def SNPfold_pool(list_snps):
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #performs the ddG calculations in parallel
    list_snps = p.map(run_SNPfold, list_snps)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()
    return list_snps

class SNP:
    def __init__(self):
        x=1

def draw_PFPs(list_cands, snp_filename, path_to_SPOT):
    for cand in list_cands:
        ct_filename = path_to_snpfold+'/'+'wt_region_'+str(cand.start)+'_'+str(cand.stop)+'.ct'
        roi_filename = '/'.join(snp_filename.split('/')[:-1])+'_ROI_'+str(cand.start)+'_'+str(cand.stop)+'.txt'
        os.system('Rscript PFP_generalized_piped_into.R %s %s %s %s' % (int(cand.start), int(cand.stop), roi_filename,ct_filename))
        os.remove(ct_filename)
        os.remove(ct_filename.replace('.ct', '.bpseq'))
        os.remove(ct_filename.replace('.ct', '.prob'))

    
    
def master(input):
    file = open(input, 'r')
    lines = file.readlines()
    list_snps = []
    header = lines[0].split('\t')
    for line in lines[1:]:
        snp = SNP()
        data = line.split('\t')
        for n, pt in enumerate(data):
            setattr(snp,header[n].replace('\n',''),pt.replace('\n','').replace('\xc2\xb1','+/-'))
        list_snps.append(snp)
    #list_snps = SNPfold_pool(list_snps)
    file.close()
    write_output(list_snps, input)
    list_cands = extract_candidates(args.c)
    analyze_candidates(list_snps,list_cands,args)
    write_candidate_outputs(list_cands, input)
    #draw_PFPs(list_cands, input)

def write_candidate_outputs(list_cands, input):
    write_master_cand_output(list_cands, input)
    write_indiv_cand_output(list_cands, input)

def write_master_cand_output(list_cands, snp_filename):
    cand_filename = '/'.join(snp_filename.split('/')[:-1])+'_Identified_PFPs.txt'
    cand_file = open(cand_filename,'w')
    cols = ['start', 'stop','label','length', 'SNP_count','rbSN_count','psSN_count', 'sequence']
    header = '\t'.join(cols)
    outstring = header +'\n'
    for cand in list_cands:
        if cand.PFP == True:
            info = cand.__dict__
            outlist = []
            for col in cols:
                outlist.append(str(info[col]))
            outline = '\t'.join(outlist) +'\n'
            if outline!='\n':
                outstring += outline
    cand_file.write(outstring)
    cand_file.close()

def write_indiv_cand_output(list_cands, snp_filename):
    for cand in list_cands:
        cand_filename = '/'.join(snp_filename.split('/')[:-1])+'_ROI_'+str(cand.start)+'_'+str(cand.stop)+'.txt'
        write_output(cand.list_snps,cand_filename)


def write_output(snps_in_candidate_regions,snp_data_libraryname):
    cols = ['pos', 'wt', 'alt','occ','mutated_protein','amino_acid','wt_aa','alt_aa','mutation_type','region','region_stop','wt_dg','alt_dg','ddg','cc_bpprob','cc_shannon','d_max','r_min','pk_bps','pk_sensitivity','pk_ppv','pb_pk_bps','pb_pk_sensitivity','pb_pk_ppv','acc']
    header = "\t".join(cols)
    outstring = header+'\n'
    for SNP in snps_in_candidate_regions:
        info = SNP.__dict__
        outlist = []
        for col in cols:
            outlist.append(str(info[col]).replace('\xb1','+/-'))
        outline = '\t'.join(outlist) +'\n'
        if outline!='\n':
            outstring += outline
    snp_data_library = open(snp_data_libraryname, 'w')
    snp_data_library.write(outstring)

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


#acceptable rate SNP: 0.3
#rate rbSN: 0.1
#rate psSN: 0.2

def analyze_candidates(list_snps_in_candidate_regions, list_candidates,args):
    apply_thresholds(list_snps_in_candidate_regions,args)
    for cand in list_candidates:
        for snp in list_snps_in_candidate_regions:
            if (str(snp.region) == str(cand.start)):
                if (snp.mutation_type == 'Missense Conservative' or snp.mutation_type == 'Silent'):
                    cand.list_snps.append(snp)
    
    for cand in list_candidates:
        cand.SNP_count = len(cand.list_snps)
        cand.rbSN_count = 0 
        cand.psSN_count = 0
        for snp in cand.list_snps:
            if (snp.S_rbSN + snp.R_rmin_rbSN + snp.R_dmax_rbSN) > 0:
                cand.rbSN_count += 1
            if snp.S_psSN + snp.P_psSN > 0:
                cand.psSN_count += 1
        if (float(cand.SNP_count/cand.length) <= 0.3) and (float(cand.rbSN_count/cand.length)<=0.1) and (float(cand.psSN_count/cand.length) <= 0.2):
            cand.PFP = True
        else:
            cand.PFP = False
        

def apply_thresholds(list_snps,args):
    psSN_threshold = float(args.t)
    snpfold_threshold = float(args.st)
    r_min_threshold = float(args.rmin)
    d_max_threshold = float(args.dmax)
    for snp in list_snps:
        snp.S_psSN = 0
        snp.P_psSN = 0
        snp.S_rbSN = 0
        snp.R_rmin_rbSN = 0
        snp.R_dmax_rbSN = 0
        if (snp.pk_sensitivity) != '':
            if float(snp.pk_sensitivity) <= psSN_threshold:
                snp.S_psSN = 1
        if (snp.pb_pk_sensitivity) != '':
            if float(snp.pb_pk_sensitivity) <= psSN_threshold:
                snp.P_psSN = 1
        if (snp.cc_bpprob) != '':
            if float(snp.cc_bpprob) <= snpfold_threshold:
                snp.S_rbSN = 1
        if (snp.r_min) != '' and (snp.r_min  != 'ERROR'):
            if float(snp.r_min) <= r_min_threshold:
                snp.R_rmin_rbSN = 1
        if (snp.d_max) != '' and (snp.d_max != 'ERROR'):
            if float(snp.d_max) <= d_max_threshold:
                snp.R_dmax_rbSN = 1

master(input)