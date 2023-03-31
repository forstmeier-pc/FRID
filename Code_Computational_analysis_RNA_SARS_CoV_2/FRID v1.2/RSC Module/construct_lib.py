import subprocess
import multiprocessing as mp
import time
import os
import argparse
import meta_pk as mpk
import sys

time_id = (time.ctime(time.time())).replace(' ', '_').replace(':', '.')

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', type=str, help='input reference sequence')
    parser.add_argument('-s', type=str, help='input RSC library search parameter')
    parser.add_argument('-d', type=str, help='input desired library depth')
    parser.add_argument('-c', type=str, help='input candidate regions in standard format')
    parser.add_argument('-g', type=str, help='gff3 file of the species (optional)')

    args = parser.parse_args()
    
    return args

def construct_library(seach, library_depth, run_ID):
    try:
        os.mkdir('RSC_sequence_library')
    except:
        pass
    os.chdir('RSC_sequence_library')
    library_fas = get_library_sequences(seach, library_depth, run_ID)
    library_drct = library_fas.replace('.fas','/')
    try:
        os.mkdir(library_drct)
    except:
        pass
    os.rename(library_fas,library_drct+library_fas)
    os.chdir(library_drct)
    extract_fasta_files(library_fas)
    os.chdir('../../')
    return 'RSC_sequence_library/'+library_drct

def extract_fasta_files(input_file):
    input_file_open = open(input_file, 'r')
    retrieved_fasta_names = open('library_fastas_QC.txt', 'a')
    input_lines = input_file_open.readlines()
    name = ' '
    sequence = ''
    for line in input_lines:
        if '>' in line and name == ' ':
            name = line
        elif '>' not in line and name != ' ':
            sequence += line
        elif '>' in line and name != ' ':
            temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','_')[:51]+'.fasta', 'w')
            temp_file.write(name+'\n'+sequence)
            temp_file.close()
            retrieved_fasta_names.write(name)
            #retrieved_fasta_names.close()
            sequence = ''
            name = line
    os.remove(input_file)
    temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','_')[:51]+'.fasta', 'w')
    temp_file.write(name+'\n'+sequence)
    temp_file.close()
    retrieved_fasta_names.write(name)

def get_library_sequences(search, library_depth, run_ID):
    subprocess.run('esearch -db nucleotide -query "%s complete" | efetch -stop %s -format fasta > library_%s.fas'%(search, str(library_depth), run_ID), shell=True)
    return 'library_%s.fas' % (run_ID)

def pool_manager_blastn(file_list):
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #performs the blastn in parallel
    p.map(run_blastn_on_sequence, file_list)
    #extracts SNPs from blastn results in parallel
    p.map(analyze_blastn_results, file_list)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()

def individual_blast(query, subject, output):
    subprocess.run('blastn -query %s -subject %s -outfmt 1 -max_hsps 1 -out %s' % (query, subject, output),shell = True)

def run_blastn_on_lib(ref_seq, path_to_RSC_library):
    file_list = []
    os.chdir(path_to_RSC_library)
    for file in os.listdir(os.getcwd()):
        if '.fasta' in file:
            file_list.append(file)
    pool_manager_blastn(file_list)
    os.chdir('../../')

def run_blastn_on_sequence(file):
    blastn_output_filename = file.replace('.fasta', '_blastn.txt')
    blastn_output_file = open(blastn_output_filename,'w')
    individual_blast(file, '../../'+ref_seq, blastn_output_filename)

def analyze_blastn_results(blastn_file):
    filename = blastn_file.replace('.fasta', '_blastn.txt')
    snp_filename = blastn_file.replace('.fasta', '_snps.txt')
    file = open(filename, 'r')
    file_lines = file.readlines()
    list_snps = []
    line_beginning_loci = 0
    #this for loop examines each line of the data file
    for number_line, line in enumerate(file_lines):
        addition = False
        line_items = line.split(' ')
        #this conditional selects for the lines of interest
        if line_items[0] == 'Subject_1':
            #this for loop goes thorguh each entry in the data line
            for num_item, item in enumerate(line_items):
                #this defines where in the genome the SNP is located
                if str(item).isdigit() == True:
                    line_beginning_loci = int(item)
                    continue
                #this finds and identifies the wt base
                if '.' in item:
                    for wt_number, wt_character in enumerate(item):
                        if wt_character in ['A', 'G', 'C', 'T']:
                            snp_loci = wt_number + line_beginning_loci
                            wt_base = wt_character
                            prev_items = file_lines[number_line-1].split(' ')
                            for item in prev_items:
                                if ('A' in item) or ('T' in item) or ('C' in item) or ('G' in item):
                                    for alt_number, alt_character in enumerate(item):
                                        if alt_number == wt_number and alt_character != '.':
                                            snp_string = wt_character+'\t'+str(snp_loci)+'\t'+alt_character
                                            list_snps.append(snp_string)
                                    break
    output_string = '\n'.join(list_snps)
    snp_file = open(snp_filename,'w')
    snp_file.write(output_string)
    snp_file.close()
    #os.remove(filename)

def collect_SNPs_to_lib(path_to_RSC_lib, run_ID):
    os.chdir(path_to_RSC_lib)
    all_lines = []
    for filename in os.listdir(os.getcwd()):
        if '_snps.txt' in filename:
            acc = filename.split('_')[0]
            file = open(filename,'r')
            lines = file.readlines()
            for line in lines:
                all_lines.append(line.replace('\n','')+'\t'+acc)
            #os.remove(filename)
    all_lines.sort(key=lambda x:int(x.split('\t')[1]))
    all_lines = count_occurrences(all_lines)
    header = '\t'.join(['wt', 'pos', 'alt', 'occ', 'acc'])
    output_string = header+'\n'+'\n'.join(all_lines)+'\n'
    output_filename = 'library_SNPs_%s.txt' % (run_ID)
    output_file = open(output_filename,'w')
    output_file.write(output_string)
    output_file.close()
    os.chdir('../../')
    return path_to_RSC_lib+output_filename

def count_occurrences(lines):
    indecies_2_remove = []
    new_lines = []
    for n, line in enumerate(lines):
        occ = 1
        info = line.split('\t')
        acc = info[3].replace('\n','')
        for n2, line2 in enumerate(lines[n+1:]):
            info2 = line2.split('\t')
            if info[:3] == info2[:3]:
                occ += 1
                acc += ', '+info2[3]
                if n+n2+1 not in indecies_2_remove:
                    indecies_2_remove.append(n+n2+1)
        new_line = '\t'.join(info[:3])+'\t'+str(occ)+'\t'+acc
        new_lines.append(new_line)
    for index in reversed(indecies_2_remove):
        del new_lines[index]
    return new_lines

class Candidate:
    def __init__(self, start, stop, sequence, label):
        self.start = int(start)
        self.stop = int(stop)
        self.label = label
        self.sequence = sequence
        self.lenth = len(sequence)

def extract_candidates(PK_candidate_file):
    region_file = open(PK_candidate_file,'r')
    lines = region_file.readlines()
    list_candidates = []
    for line in lines[1:]:
        info = line.split('\t')
        region = Candidate(info[0], info[1], info[2], info[3])
        list_candidates.append(region)
    return list_candidates

class SNP:
    def __init__(self, wt, pos, alt, occ, acc):
        self.wt = wt
        self.alt = alt
        self.pos = int(pos)
        self.occ = int(occ)
        self.acc = acc.replace('\n','')
        self.ID = wt+str(pos)+alt

def extract_snps(path_to_snp_lib):
    snp_lib_file = open(path_to_snp_lib,'r')
    lines = snp_lib_file.readlines()
    list_snps = []
    for line in lines[1:]:
        info = line.split('\t')
        snp = SNP(info[0],info[1], info[2], info[3], info[4])
        list_snps.append(snp)
    return list_snps

def get_snps_in_candidate_regions(PK_candidate_file, path_to_snp_lib):
    list_candidates = extract_candidates(PK_candidate_file)
    list_snps = extract_snps(path_to_snp_lib)
    list_cand_snps = []
    for cand in list_candidates:
        for snp in list_snps:
            if snp.pos < cand.stop and snp.pos > cand.start:
                if snp.alt.capitalize() in ['A','C','T','G']:
                    snp.region = cand.start
                    list_cand_snps.append(snp)
    return list_cand_snps

def make_wt_and_alt_sequence(SNP):
    wt_filename = 'wt_'+SNP.ID+'.fasta'
    wt_file = open(wt_filename,'w')
    alt_filename = 'alt_'+SNP.ID+'.fasta'
    alt_file = open(alt_filename, 'w')
    if ref_seq_string[SNP.pos-1] == SNP.wt:
        five_prime_flanker = -61
        if SNP.pos <= 61:
            five_prime_flanker = SNP.pos-1
        ref_seq_window = ref_seq_string[int(SNP.pos+five_prime_flanker):SNP.pos+60]
        wt_file.write('>wt_'+SNP.ID+'\n'+ref_seq_window)
        alt_seq_window = ref_seq_string[int(SNP.pos)+five_prime_flanker:int(SNP.pos)-1] + SNP.alt+ref_seq_string[int(SNP.pos):int(SNP.pos)+60]
        alt_file.write('>alt_'+SNP.ID+'\n'+alt_seq_window)
        SNP.wt_window = ref_seq_window
        SNP.alt_window = alt_seq_window
        return [wt_filename, alt_filename]
    else:
        print('WT BASE DOES NOT MATCH SNP BASE')

def get_individual_ddg(SNP):
    filenames = make_wt_and_alt_sequence(SNP)
    wt_filename = filenames[0]
    alt_filename = filenames[1]
    SNP.wt_dg = run_RNAstructure_fold(wt_filename, SNP.ID)
    SNP.alt_dg = (run_RNAstructure_fold(alt_filename, SNP.ID))
    alt_fe = float(SNP.alt_dg.split(' ')[0])
    alt_fe_er = float(SNP.alt_dg.split(' ')[2])
    wt_fe = float(SNP.wt_dg.split(' ')[0])
    wt_fe_er = float(SNP.wt_dg.split(' ')[2])
    total_err = round(wt_fe_er + alt_fe_er, 2)
    d_fe = round(alt_fe - wt_fe, 2)
    SNP.ddg = (str(d_fe)+' '+SNP.wt_dg.split(' ')[1]+' '+str(total_err))
    os.remove(wt_filename)
    os.remove(alt_filename)
    return SNP

def run_RNAstructure_fold(seq_filename, ID):
    fold_output_name = 'Fold_output_%s.ct' % (ID)
    fold_output = open(fold_output_name, 'w')
    subprocess.run('Fold -mfe %s %s' % (seq_filename, fold_output_name),shell=True)
    efn2_output_name = 'efn2_output_%s.txt' % (ID)
    efn2_output = open(efn2_output_name, 'w')
    subprocess.run('efn2 %s %s' % (fold_output_name, efn2_output_name), shell=True)
    efn2_reader = open(efn2_output_name, 'r', encoding='utf-8')
    mfe_lines = efn2_reader.readlines()
    info_string = mfe_lines[0].split('= ')[1].replace('\n','')
    os.remove(fold_output_name)
    os.remove(efn2_output_name)
    return info_string    

def get_ddgs(snps_in_candidate_regions, ref_seq_string):
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #performs the ddG calculations in parallel
    snps_in_candidate_regions = p.map(get_individual_ddg, snps_in_candidate_regions)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()
    return snps_in_candidate_regions

def get_RNAsnp(snps_in_candidate_regions):
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #p = mp.Pool(1)
    #performs the ddG calculations in parallel
    snps_in_candidate_regions = p.map(run_RNAsnp, snps_in_candidate_regions)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()
    return snps_in_candidate_regions

def make_snp_file(SNP):
    snp_filename = SNP.ID +'.snp'
    snp_file = open(snp_filename, 'w')
    snp_file.write(SNP.wt+str(61)+SNP.alt)
    snp_file.close()
    return snp_filename

def run_RNAsnp(SNP):
    filenames = make_wt_and_alt_sequence(SNP)
    wt_filename = filenames[0]
    alt_filename = filenames[1]

    RNAsnp_output_name = 'RNAsnp_output_%s.txt' % (SNP.ID)
    RNAsnp_output = open(RNAsnp_output_name, 'w')
    
    SNP_file = make_snp_file(SNP)
    subprocess.run('RNAsnp -f %s -s %s -w 100 > %s' % (wt_filename, SNP_file, RNAsnp_output_name), shell=True)
    d_max = False
    r_min = False
    change = False
    RNAsnp_output.close()
    RNAsnp_output = open(RNAsnp_output_name, 'r')
    RNAsnp_output_lines = RNAsnp_output.readlines()
    for num_line, line in enumerate(RNAsnp_output_lines):
        split_lines = line.split('\t')
        if 'SNP' in split_lines[0]:
            try:
                data_line = RNAsnp_output_lines[num_line+1].split('\t')
                break
            except:
                error_RNAsnp_file = open('RNAsnp_error_log.txt', 'a')
                error_RNAsnp_file.write(RNAsnp_output_name+'\n')
                os.remove(RNAsnp_output_name)
                return SNP
    if len(data_line) == 1:
        error_RNAsnp_file = open('RNAsnp_error_log.txt', 'a')
        error_RNAsnp_file.write(RNAsnp_output_name+'\n')
        os.remove(RNAsnp_output_name)
        return SNP
    
    SNP.r_min = data_line[6]
    SNP.d_max = data_line[9]
    os.remove(RNAsnp_output_name)
    return SNP

def get_pseudoSNitches(threshold, snps_in_candidate_regions,path_to_SPOT_RNA):
    for SNP in snps_in_candidate_regions:
        os.chdir(path_to_SPOT_RNA)
        filenames = make_wt_and_alt_sequence(SNP)
        wt_filename = filenames[0]
        alt_filename = filenames[1]
        spot_rna_output_path_huh= ''
        subprocess.run('python3 SPOT-RNA.py --inputs '+str(wt_filename)+' --outputs ./', shell=True)
        subprocess.run('python3 SPOT-RNA.py --inputs '+str(alt_filename)+' --outputs ./', shell=True)
        comparison = mpk.scorer_for_pks(SNP.region, 'alt_'+SNP.ID+'.ct', 'wt_'+SNP.ID+'.ct', 2,2)
        os.remove(wt_filename)
        os.remove(alt_filename)
        os.remove(wt_filename.replace('.fasta', '.bpseq'))
        os.remove(wt_filename.replace('.fasta', '.ct'))
        os.remove(wt_filename.replace('.fasta', '.prob'))
        os.remove(alt_filename.replace('.fasta', '.bpseq'))
        os.remove(alt_filename.replace('.fasta', '.ct'))
        os.remove(alt_filename.replace('.fasta', '.prob'))
        os.chdir('~/scratch')
        SNP.pk_bps = comparison[5]
        SNP.pk_sensitivity = comparison[4]/comparison[5]
        SNP.pk_ppv = comparison[6]/comparison[6]
    return snps_in_candidate_regions

def set_up_intermediate_for_SNPfold(snps_in_candidate_regions,RSC_library_path, run_ID):
    snp_data_libraryname = RSC_library_path +'library_SNPS_data_'+run_ID+'.txt'
    snp_data_library = open(snp_data_libraryname, 'w')
    cols = []
    for col in snps_in_candidate_regions[0].__dict__.keys():
        cols.append(col)
    header = "\t".join(cols)
    outstring = header+'\n'
    for SNP in snps_in_candidate_regions:
        info = SNP.__dict__
        outlist = []
        for col in cols:
            outlist.append(str(info[col]))
        outline = '\t'.join(outlist) +'\n'
        outstring += outline
    snp_data_library.write(outstring)
    return snp_data_libraryname

def get_ref_seq_string(ref_seq):
    ref_file = open(ref_seq,'r')
    lines = ref_file.readlines()
    seq_string = ''
    for line in lines[1:]:
        seq_string += line.replace('\n','')
    return seq_string

class GFF3_LINE:
    def __init__(self, subset, start, end, name):
        self.subset = subset
        self.start = int(start)
        self.end = int(end)
        self.name = name

def get_gff3_line(line):
    line = line.split('\t')
    notes = line[8].split(';')
    product = 'NULL'
    for note in notes:
        if 'product' in note:
            product = note.split('product=')[1]
    return GFF3_LINE(line[2],line[3],line[4], product)

def get_gff_lines(gff_filename):
    gff3_file = open(gff_filename,'r')
    gff3_lines = gff3_file.readlines()
    gff3_lines_data = []
    for line in gff3_lines[7:-1]:
        gff3_lines_data.append(get_gff3_line(line))
    return gff3_lines_data
    
    #use this to extract lines of the gff file to spread over the pool of labeling threads

def label_AAs(snps_in_candidate_regions,gff_filename):
    if gff_filename == None:
        return snps_in_candidate_regions
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #performs the ddG calculations in parallel
    snps_in_candidate_regions = p.map(label_AA_of_SNP, snps_in_candidate_regions)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()
    return snps_in_candidate_regions

def fetch_sequence(start, end_inclusive):
    desired_sequence = ref_seq_string[start-1:end_inclusive]
    return desired_sequence

def label_AA_of_SNP(SNP):
    codon_filename = 'codon_table.txt'
    codon_file = open(codon_filename,'r')

    codon_lines = codon_file.readlines()
    codon_dict = {}
    codon_dict_single_letter = {}
    for line in codon_lines[2:]:
        data = line.split('\t')
        codon_dict[data[0]] = data[1]
        codon_dict_single_letter[data[1]] = data[2]
    red_aa = ['A','V','F','P','M','I','L','W'] #small hydrophobic
    blue_aa = ['D','E'] #acidic
    magenta_aa = ['R','K'] # Basic
    green_aa = ['S','T','Y','H','C','N','G','Q'] #hydroxyl, sulfhydryl,amine, G
    list_colors_aas = [red_aa,blue_aa,magenta_aa,green_aa]

    mutation = ''
    mutation_type = ''
    mutated_protein = ''
    protein_mutation_location_in_protein = ''
    output_line = line
    wt_aa = ''
    mu_aa = ''
    for gff in gff_obj_list:
        if gff.start < int(SNP.pos) < gff.end:
            if 'CDS' in gff.subset:
                protein_mutation_location_in_protein = ((int(SNP.pos)-gff.start)//3)+1
                
                if ((int(SNP.pos)-gff.start+1)%3==0):
                    codon = fetch_sequence(int(SNP.pos)-2,int(SNP.pos))
                    mu_codon = codon[:2]+SNP.alt.capitalize()
                    wt_aa =codon_dict[codon]
                    if 'Del' in mu_codon:
                        mu_aa = 'Frameshift'
                    else:
                        mu_aa = codon_dict[mu_codon]

                if ((int(SNP.pos)-gff.start+1)%3==2):
                    codon = fetch_sequence(int(SNP.pos)-1,int(SNP.pos)+1)
                    mu_codon = codon[0]+SNP.alt.capitalize()+codon[2]
                    wt_aa =codon_dict[codon]
                    if 'Del' in mu_codon:
                        mu_aa = 'Frameshift'
                    else:
                        mu_aa = codon_dict[mu_codon]

                if ((int(SNP.pos)-gff.start+1)%3==1):
                    codon = fetch_sequence(int(SNP.pos),int(SNP.pos)+2)
                    wt_aa =codon_dict[codon]
                    mu_codon = SNP.alt.capitalize()+codon[1:]
                    if 'Del' in mu_codon:
                        mu_aa = 'Frameshift'
                    else:
                        mu_aa = codon_dict[mu_codon]

                if mu_aa != 'Frameshift':
                    mutation = '%s>%s' % (codon_dict_single_letter[wt_aa], codon_dict_single_letter[mu_aa])
                else:
                    mutation = '%s>%s' % (codon_dict_single_letter[wt_aa], 'FS')
                mutation_type = ''
                mutated_protein = gff.name
                protein_mutation_location_in_protein = str(protein_mutation_location_in_protein)

                if mu_aa == 'Frameshift':
                    mutation_type = 'Frameshift'
                elif codon_dict_single_letter[mu_aa] == 'O':
                    mutation_type ='Nonsense'
                elif mu_aa != wt_aa:
                    for lst in list_colors_aas:
                        if (codon_dict_single_letter[wt_aa] in lst) and (codon_dict_single_letter[mu_aa] not in lst):
                            mutation_type='Missense Non-conservative'
                            break
                        elif (codon_dict_single_letter[wt_aa] in lst) and (codon_dict_single_letter[mu_aa] in lst):
                            mutation_type ='Missense Conservative'
                            break

                else:
                    mutation_type='Silent'
                break

    SNP.mutated_protein = mutated_protein
    SNP.mutation_type = mutation_type
    SNP.amino_acid = protein_mutation_location_in_protein
    SNP.wt_aa = wt_aa
    SNP.alt_aa = mu_aa
    return SNP

def master():
    threshold = 0
    args = get_args()
    search = args.s
    library_depth = args.d
    candidate_filename = args.c
    gff_filename = args.g
    global ref_seq
    ref_seq = args.r
    global ref_seq_string
    ref_seq_string = get_ref_seq_string(ref_seq)
    run_ID = search.replace(' ','_')+'_'+str(library_depth)+'_'+time_id
    RSC_library_path = construct_library(search,library_depth,run_ID)
    run_blastn_on_lib(ref_seq, RSC_library_path)
    SNP_lib_path = collect_SNPs_to_lib(RSC_library_path, run_ID)
    snps_in_candidate_regions = get_snps_in_candidate_regions(candidate_filename, SNP_lib_path)
    snps_in_candidate_regions = get_ddgs(snps_in_candidate_regions, ref_seq_string)
    snps_in_candidate_regions = get_pseudoSNitches(threshold, snps_in_candidate_regions, '../programs/SPOT-RNA')
    snps_in_candidate_regions = get_RNAsnp(snps_in_candidate_regions)
    global gff_obj_list
    gff_obj_list = get_gff_lines(gff_filename)
    #snps_in_candidate_regions = label_AAs(snps_in_candidate_regions, gff_filename)
    intermediate_filename = set_up_intermediate_for_SNPfold(snps_in_candidate_regions,RSC_library_path, run_ID)
    return intermediate_filename

file = master()
sys.exit(file)

def run_blastn_orig(PK_candidates_input_file, ref_seq):
   
    #this is the working directory of the pipeline
    os.chdir('SARS_CoV_2')
    #The ref seq of SARS-CoV-2
    ref_file = ref_seq
    #this directory contains isolate genomes of the virus
    os.chdir('all_SARS_CoV_2_genomes_no_ref')
    all_snps = []
    all_files = []
    #for file in os.listdir():
    #    if 'Pass' not in file and 'blastn' not in file and 'snps' not in file:
    #        extract_fasta_files(file)
    
    #this loops through every genome in the file
    for num_file, file in enumerate(os.listdir()):
        print('\n==========================================================\n'+file+'\n')
        #this ensures that it is a genomic file that is being scanned
        if 'temp_file' not in file and 'Pass' not in file and 'blastn' not in file and 'snps' not in file and 'priority' in file:
        #if file == 'MT451283.1Severeacuterespiratorysyndromecoronavirus.fasta':
            all_files.append(file)
            snp_file_file_name = file.replace('.fasta', '_snps.txt')
            snp_file_file = open(snp_file_file_name, 'w')
            snp_file_string = ''
            print(file)
            temp_file_name = file.replace('.fasta','_blastn.txt')
            temp_file = open(temp_file_name, 'w')
            temp_file.close()
            #this blasts the sequence
            #individual_blast(file, ref_file, temp_file_name)
            subprocess.run('blastn -query %s -subject %s -outfmt 1 -out %s' % (file, ref_file, temp_file_name),shell = True)
            os.chdir('../')
            ###directory location is in SARS_CoV_2
            #this runs the rest of the script
            file_snps = check_snps_in_pks('all_SARS_CoV_2_genomes_no_ref/'+temp_file_name, putative_pks_input_file, this_run_num)
            #this collects the sequence specific SNPs
            file_snps.sort(key=lambda x:x.pos)
            for lst in file_snps:
                all_snps.append(lst)
                snp_file_string += (str(lst.wt) +' '+str(lst.pos)+' '+str(lst.alt)+'\n')
            snp_file_file.write(file+'\n'+snp_file_string)
            print(os.getcwd(), 'CCCC')
            os.chdir(home_directory_path+'SARS_CoV_2/all_SARS_CoV_2_genomes_no_ref')
        else:
            continue
    for num_file, file in enumerate(os.listdir()):
        if num_file == 200:
            break
        print('\n==========================================================\n'+file+'\n')
        #this ensures that it is a genomic file that is being scanned
        if 'temp_file' not in file and 'Pass' not in file and 'blastn' not in file and 'snps' not in file:
        #if file == 'MT451283.1Severeacuterespiratorysyndromecoronavirus.fasta':
            all_files.append(file)
            snp_file_file_name = file.replace('.fasta', '_snps.txt')
            snp_file_file = open(snp_file_file_name, 'w')
            snp_file_string = ''
            print(file)
            temp_file_name = file.replace('.fasta','_blastn.txt')
            temp_file = open(temp_file_name, 'w')
            temp_file.close()
            #this blasts the sequence
            #individual_blast(file, ref_file, temp_file_name)
            subprocess.run('blastn -query %s -subject %s -outfmt 1 -out %s' % (file, ref_file, temp_file_name),shell = True)
            os.chdir('../')
            ###directory location is in SARS_CoV_2
            #this runs the rest of the script
            file_snps = check_snps_in_pks('all_SARS_CoV_2_genomes_no_ref/'+temp_file_name, putative_pks_input_file, this_run_num)
            #this collects the sequence specific SNPs
            file_snps.sort(key=lambda x:x.pos)
            for lst in file_snps:
                all_snps.append(lst)
                snp_file_string += (str(lst.wt) +' '+str(lst.pos)+' '+str(lst.alt)+'\n')
            snp_file_file.write(file+'\n'+snp_file_string)
            print(os.getcwd(), 'CCCC')
            os.chdir(home_directory_path+'SARS_CoV_2/all_SARS_CoV_2_genomes_no_ref')
        else:
            continue
    store_rando_files(all_files)

    #this collects all the snps from the different sequences and makes sure there are no repeats and counts the number of specific repeats there are
    tracking_snps = []
    for place, snp in enumerate(all_snps):
        check_snp = snp.wt+str(snp.pos)+snp.alt
        if check_snp not in tracking_snps:
            tracking_snps.append(check_snp)
            occurence = 1
            for other_place, other_snps in enumerate(all_snps):
                other_check = other_snps.wt +str(other_snps.pos)+other_snps.alt
                if check_snp == other_check and place != other_place:
                    occurence += 1
                    snp.isolate.append(other_snps.isolate[0])
            snp.occurence = occurence
    snps_for_deletion = []
    for place, snp in enumerate(all_snps):
        if snp.occurence == 0:
            snps_for_deletion.append(place)

    i = 1
    while i <= len(snps_for_deletion):
        del all_snps[snps_for_deletion[-i]]
        i += 1


    time_id = (time.ctime(time.time())).replace(' ', '_').replace(':', '.')
    all_snps_file = open('all_snps_%s_%s.txt' % (str(this_run_num), time_id), 'w')
    all_snps_string = 'loci\twt\tmu\tOccurences\tdG SNP\tdG wt\tddG\trbSN\tpsSN\tisolates\n'
    all_snps_name_string = ''
    all_snps.sort(key=lambda x:x.pos)
    intermediate_snpfold_storage_file = open('intermediate.txt', 'a')
    for snp in all_snps:
        snp.initialize_list()
        new_snp = [snp.wt, snp.alt, str(snp.pos),str(snp.pk)]
        if snp.pk !=0:
            intermediate_snpfold_storage_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n '%(str(snp.pk), str(snp.snp_string), snp.ref_pk_filename, snp.file, str(snp.rnasnp_results), this_run_num, str(new_snp),'riboSNitch_programs_data', 'True'))
      
        #try:
        print(snp.extol(), 'final_check')
        print(snp.isolate)
        isolates = ', '.join(snp.isolate)
        all_snps_string += ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (str(snp.pos), snp.wt, snp.alt,str(snp.occurence), snp.snp_dg, snp.wt_dg, snp.ddg, snp.rnasnp_results, snp.pk_results, isolates))
            #all_snps_string += (snp[1]+'\t'+str(snp[2])+'\t'+snp[3]+ '\tOccurences: '+str(snp[10])+'\tdG SNP: '+str(snp[7])+'\tdG WT: '+str(snp[8])+'\tddG: '+str(snp[9])+'\trbSN: '+str(snp[5])+'\tpsSN: '+str(snp[6])+'\n')
        #except:
        #    print("BAD_WRITING_BAD_WRITING_ADMIN")
            #all_snps_string += (snp[1]+'\t'+str(snp[2])+'\t'+snp[3]+'\tOccurences: '+str(snp[8])+'\tdG SNP: '+str(snp[5])+'\tdG WT: '+str(snp[6])+'\tddG: '+str(snp[7])+'\n')
    intermediate_snpfold_storage_file.close()
    for name in all_files:
        all_snps_name_string += name+'\n'
    all_snps_file.write(time_id+'\n'+all_snps_name_string+'\n'+all_snps_string)
    all_snps_file.close()
    organize_genome_folder_SARS()
    #list_pks = extract_putative_pks(home_directory_path+'SARS_CoV_2/'+ putative_pks_input_file)
    #compile_ROI_data(this_run_num, list_pks)
    #compare_ROIs_to_control(list_pks)
    os.rename(home_directory_path+'SARS_CoV_2/all_SARS_CoV_2_genomes_no_ref/intermediate.txt', home_directory_path+'SARS_CoV_2/intermediate.txt')
    print(time.ctime(time.time()), 'time')

#this function takes the extracted data from both the pks and the blast to find SNPs in pks
def check_snps_in_pks(input_file, putative_pks_input_file, this_run_num):
    
    
    #this grabs all the SNPs found by the BLASTn search
    list_snps = extract_data_from_blastn(input_file)
    #this grabs the regions of interest from the input file
    list_pks = extract_putative_pks(putative_pks_input_file)
    try:
        os.mkdir('putative_SNitches')
    except:
        pass
    #looks through every pk
    for pks in list_pks:
        SNPS_in_pk = []
        indel_questino = False
        #for every pk, this loop looks through all SNPs to see if any fall into this region
        for num_SNP, SNP in enumerate(list_snps):
            if SNP.wt == SNP.alt:
                list_snps.remove(SNP)
        for SNP in list_snps:
            if SNP.pos in range(int(pks.five_prime_boundary), int(pks.three_prime_boundary)):
                #this conditional checks whether an indel is also in this region
                if 'DEL' in SNP.wt or '-' in SNP.wt:
                    indel_questino = True
                    try:
                        indel(input_file, SNP)
                    except:
                        pass
                    continue
                SNPS_in_pk.append(SNP)
        if indel_questino == True:
            break
        #loops through only SNPs in regions of interest
        for SNP in SNPS_in_pk:
            print(os.getcwd())
            os.chdir(home_directory_path+'SARS_CoV_2/putative_SNitches')
            
            #now that SNPs have been filtered, this runs all the 
            print(SNP)
            snp_pks = pks.create_SNP(SNP)
            if snp_pks == False:
                snp_construction_error_file = open('snp_construction_error_file', 'a')
                snp_construction_error_file.write(input_file+' '+str(SNP.extol())+'\n')
                continue
            snp_pks.snp_string = SNP.wt+str(SNP.pos-int(pks.five_prime_boundary)+1)+SNP.alt
            SNP.snp_string = snp_pks.snp_string
            SNP.pk = snp_pks.five_prime_boundary
            degree_support = 0

            print ('\nAnalyzing SNP with SPOT-RNA...')
            if snp_pks.snp_string[0] == '-':
                continue
            spot_score = run_spot_snp(snp_pks.five_prime_boundary, snp_pks.five_bond, snp_pks.snp_string, snp_pks.sequence, degree_support, input_file)
            #we are in putativeSNitches 
            loci = SNP.pos-int(pks.five_prime_boundary)
            snp_loci = 61
            if loci < 61:
                snp_pks.relevant_sequence = snp_pks.sequence[:loci+60]
                pks.relevant_sequence = pks.sequence[:loci+60]
                snp_loci = loci
            elif loci + 60 > len(pks.sequence):
                snp_pks.relevant_sequence = snp_pks.sequence[loci-61:]
                pks.relevant_sequence = pks.sequence[loci-61:]
                snp_loci = 61
            else:
                snp_pks.relevant_sequence = snp_pks.sequence[loci-61:loci+60]
                pks.relevant_sequence = pks.sequence[loci-61:loci+60]
                snp_loci = 61
            print(loci, SNP.extol(), 'OOOO')
            print(snp_pks.sequence)
            print(snp_pks.relevant_sequence, 'NNNNN')
            print(pks.sequence)
            print(pks.relevant_sequence, 'MMMMM')
            snp_pks.snp_string = SNP.wt +str(snp_loci+1) + SNP.alt
            SNP.snp_string = snp_pks.snp_string
            #try:
            #    print('Base Pair Sensitivity: %d / %d = %d \nBase Pair PPV: %d / %d = %d \nPseudoknot Sensitivity: %d / %d = %d \nPseudoknot PPV: %d / %d = %d \n' % (spot_score[0], spot_score[1], float(spot_score[0]/spot_score[1]), spot_score[2], spot_score[3], float(spot_score[2]/spot_score[3]), spot_score[4], spot_score[5], float(spot_score[4]/spot_score[5]),spot_score[6], spot_score[7], float(spot_score[6]/spot_score[7])))
            #except:
            #    print('Base Pair Sensitivity: %d / %d\nBase Pair PPV: %d / %d\nPseudoknot Sensitivity: %d / %d\nPseudoknot PPV: %d / %d\n' % (spot_score[0], spot_score[1], spot_score[2], spot_score[3], spot_score[4], spot_score[5], spot_score[6], spot_score[7]))
            rnasnp_results = riboSNitch_programs(pks, snp_pks, SNP, degree_support, input_file, this_run_num)
            print('Analyzing wt and SNP with ProbKnot...')
            prob_results = run_probknot(pks, snp_pks, SNP, input_file)

            riboSN_results = int(rnasnp_results)
            
            SNP.rnasnp_results = riboSN_results

            num_pk_prgs = 0
            pk_prgs = ''
            if spot_score != None:
                num_pk_prgs += 1
                pk_prgs = spot_score
            if prob_results != None:
                num_pk_prgs += 1
                pk_prgs = prob_results
            if spot_score != prob_results and num_pk_prgs == 2:
                pk_prgs = 'change'
            pk_string = str(num_pk_prgs)+' '+str(pk_prgs)

            SNP.pk_results = pk_string

            #run_ipknot(pks, snp_pks, SNP, input_file)
            print('\n_____________________________________\n')
            remove_excess()
            os.chdir('../')
            #Now again in SARS_CoV_2 drct
        #this conditional runs the same program if there are 2 SNPs in the same region
        #if len(SNPS_in_pk) == 2:
            #os.chdir('putative_SNitches')
            #orig_location = []
            #for snp_check in SNPS_in_pk:
            #    orig_location.append(snp_check.pos)
            #if int(orig_location[0]) - int(orig_location[1]) == 0:
            #    continue
            #snp_string_pre = ''
            #snp_pk = pks
            #for SNP in SNPS_in_pk:
            #    snp_pk = snp_pk.create_SNP(SNP)
            #    if snp_pks == False:
            #        snp_construction_error_file = open('snp_construction_error_file', 'a')
            #        snp_construction_error_file.write(input_file+' '+SNP.extol()+'\n')
            #        continue
            #adjustmenet = int(list_snps[0].pos)-61
            #snp_string_pre += list_snps[0][1]+str(61)+list_snps[0][3]+'_'+list_snps[1][1]+str(int(list_snps[1][2])-adjustmenet)+list_snps[1][3]
            #snp_pk.snp_string = snp_string_pre[:-1]
            #snp_pk.is_snp = True
            #print(snp_string_pre)
            #print(snp_pk.sequence)
            #degree_support = 0
            #print ('\nAnalyzing SNP with SPOT-RNA...')
            #spot_score = run_spot_snp(snp_pk.five_prime_boundary, snp_pk.five_bond, snp_pk.snp_string, snp_pk.sequence, degree_support, input_file)
            #loci = SNP[2]-int(pks.five_prime_boundary)
            #if loci < 61:
            #    snp_pks.sequence = snp_pks.sequence[:loci+60]
            #    pks.relevant_sequence = pks.sequence[:loci+60]
            #elif loci + 60 > len(pks.sequence):
            #    snp_pks.sequence = snp_pks.sequence[loci-61:]
            #    pks.relevant_sequence = pks.sequence[loci-61:]
            #else:
            #    snp_pks.sequence = snp_pks.sequence[loci-61:loci+60]
            #    pks.relevant_sequence = pks.sequence[loci-61:loci+60]            
            #try:
            #    print('Base Pair Sensitivity: %d / %d = %d \nBase Pair PPV: %d / %d = %d \nPseudoknot Sensitivity: %d / %d = %d \nPseudoknot PPV: %d / %d = %d \n' % (spot_score[0], spot_score[1], float(spot_score[0]/spot_score[1]), spot_score[2], spot_score[3], float(spot_score[2]/spot_score[3]), spot_score[4], spot_score[5], float(spot_score[4]/spot_score[5]),spot_score[6], spot_score[7], float(spot_score[6]/spot_score[7])))
            #except:
            #    print('Base Pair Sensitivity: %d / %d\nBase Pair PPV: %d / %d\nPseudoknot Sensitivity: %d / %d\nPseudoknot PPV: %d / %d' % (spot_score[0], spot_score[1], spot_score[2], spot_score[3], spot_score[4], spot_score[5], spot_score[6], spot_score[7]))
            #riboSNitch_programs(pks, snp_pk, SNP, degree_support, input_file, this_run_num)
            #print('Analyzing wt and SNP with ProbKnot...')
            #prob_results = run_probknot(pks, snp_pk, SNP, input_file)
            #run_ipknopks, snp_pks, SNP, input_file)

            #print('\n_____________________________________\n')
            ##PATH
            #os.chdir('../')

            
        if len(SNPS_in_pk) >= 2:
            print('\n\n\n &&&&&&&&&&&&&&&&&&&&&&          ADMIN WARNING: MANY SNPs in %s putative pseudoknot          &&&&&&&&&&&&&&&&&&&&&&n\n\n\n'%(pks.five_prime_boundary))
            multiple_snps_in_pk_region(input_file, SNPS_in_pk)
    for SNP in list_snps:
        skip = False
        if 'DEL' in SNP.wt or 'del' in SNP.alt:
            continue
        found = False
        if SNP.wt == '-':
            continue
        #for pk in list_pks:
         #   try:
          #      if int(SNP[2])-int(pk.five_prime_boundary)+1 <= 0:
          #          continue
         #       snp_string = SNP[1]+str(int(SNP[2])-int(pk.five_prime_boundary)+1)+SNP[3]
        #        snp_test_file = open('putative_SNitches/'+str(pk.five_prime_boundary)+snp_string+'_SNP.fasta', 'r')
        #        snp_test_file.close()
        #        print(str(pk.five_prime_boundary), 'FFFF')
        #        SNP.append(run_RNAstructure_fold('putative_SNitches/'+str(pk.five_prime_boundary)+snp_string+'_SNP.fasta', SNP))
        #        temp_file = open('ref.txt', 'w')
        #        temp_file.write('>ref\n'+pk.sequence[SNP[2]-61:SNP[2]+60])
        #        temp_file.close()
        #        print('here')
        #        SNP.append(run_RNAstructure_fold(temp_file, SNP))
        #        os.remove(temp_file)
        #        found = True
        #        break
        #    except:
        #        pass
        #if found == False:
        temp_name = sequence_remaining_snps(SNP)
        print(temp_name, 'AAAAA')
        SNP.snp_dg = (run_RNAstructure_fold(temp_name, SNP))
        temp_name = sequence_remaining_wt_ref(SNP)
        print(temp_name, 'BBBBB')
        SNP.wt_dg = (run_RNAstructure_fold(temp_name, SNP))
        os.remove(temp_name)
        snp_fe = float(SNP.snp_dg.split(' ')[0])
        snp_fe_er = float(SNP.snp_dg.split(' ')[2])
        wt_fe = float(SNP.wt_dg.split(' ')[0])
        wt_fe_er = float(SNP.wt_dg.split(' ')[2])
        total_err = round(wt_fe_er + snp_fe_er, 2)
        d_fe = round(snp_fe - wt_fe, 2)
        SNP.ddg = (str(d_fe)+' '+SNP.wt_dg.split(' ')[1]+' '+str(total_err))

    return list_snps



#This function extracts and formats the data from blastn, specificallly "Query-anchored showing identities" format
def extract_data_from_blastn(input_file):
    file = open(input_file, 'r')
    file_lines = file.readlines()
    list_location_snp = []
    temp_location_snp = []
    line_beginning_loci = 0
    #this for loop examines each line of the data file
    for number_line, line in enumerate(file_lines):
        addition = False
        line_items = line.split(' ')
        #this conditional selects for the lines of interest
        if line_items[0] == 'Subject_1':
            #this for loop goes thorguh each entry in the data line
            for num_item, item in enumerate(line_items):
                #this defines where in the genome the SNP is located
                if str(item).isdigit() == True:
                    line_beginning_loci = int(item)
                    continue
                #this finds and identifies the wt base
                if '.' in item:
                    for number, character in enumerate(item):
                        if character in ['A', 'G', 'C', 'T', '-']:
                                snp_loci = number + line_beginning_loci
                                temp_location_snp = [number, character, snp_loci, line_beginning_loci]
                                prev_items = file_lines[number_line-1].split(' ')
                                for item in prev_items:
                                    if ('A' in item) or ('T' in item) or ('C' in item) or ('G' in item):
                                        #for snpoly in temp_location_snp:
                                        snpoly = temp_location_snp
                                        for number, character in enumerate(item):
                                            if number == snpoly[0] and character != '.':
                                                if snp_already_analyzed(input_file,snpoly[1], snpoly[2], character) == False:
                                                    temp_snp = SNP(snpoly[1], character, snpoly[2], input_file.replace('_blastn.txt',''), snpoly[3], snpoly[0])
                                                    list_location_snp.append(temp_snp)
                                                    #list_location_snp.append([snpoly[0], snpoly[1], snpoly[2], character, snpoly[3]])
                                        break
 
                    #this loop detects insertions of the subject sequence
                    for num_item2, item2 in enumerate(file_lines[number_line+1].split(' ')):
                        if '\\' in item2:
                            lstr = file_lines[number_line+3].split(' ')
                            #' '.join(lstr).split()
                            for num_item3, item3 in enumerate(lstr):
                                item3 = item3.replace('\n', '')
                                if item3 in ['A', 'C', 'G', 'T']:
                                    deletion = True
                                    position_tracker = -1
                                    please_break  = False
                                    for final_item in line_items:
                                        if please_break == True:
                                            break
                                        if '.' not in final_item:
                                            position_tracker += 1
                                        for i in final_item:
                                            if i == '.':
                                                please_break = True
                                                break
                                            position_tracker += 1
                                    location = num_item3-position_tracker
                                    if [location, 'DEL '+item3, line_beginning_loci+location, 'del',line_beginning_loci] not in list_location_snp:
                                        if snp_already_analyzed(input_file, 'DEL '+item3,line_beginning_loci+location, character) == False:
                                            temp_snp = SNP( 'DEL '+item3, 'del', line_beginning_loci+location,  input_file.replace('_blastn.txt',''), line_beginning_loci, location)
                                            list_location_snp.append(temp_snp)
                                            
    add_snps_nomenclature(list_location_snp)
    return list_location_snp


def add_snps_nomenclature(list_snps):
    list_names = ['R', 'Y', 'M', 'K', 'S','W','H','B','V','D','N']
    list_replacements = [['G', 'A'], ['T', 'C'], ['A', 'C'], ['G', 'T'], ['G', 'C'],['A', 'T'],['A','C', 'T'], ['G', 'T', 'C'], ['G','C', 'A'],['G','T','A'], ['G','C','A','T']]
    i = 0
    #loops through the list of SNPs to determine if they contain ambiguity
    while i <= len(list_snps)-1:
        list_snps[i].alt = list_snps[i].alt.capitalize()
        #finds ambiguous bases
        if list_snps[i].alt.capitalize() in list_names:
            index_value = list_names.index(list_snps[i].alt.capitalize())
            #grabs bases with unknown identities
            if list_snps[i].alt.capitalize() == 'N':
                #print(list_replacements[-1], 'a')
                sub_list = list_replacements[-1]
                #print(sub_list, 'b')
                sub2 = sub_list
                #print(list_snps[i].wt, list_snps[i].alt, str(list_snps[i].pos))
                if list_snps[i].wt != '-':
                    sub2.remove(list_snps[i].wt)
                #print(sub2, 'c')
                list_snps[i].alt = sub2[0]
                for replacement in sub2[1:]:
                    temp_snp = SNP(list_snps[i].wt, replacement,list_snps[i].pos,'bogus/bogus',list_snps[i].line_beginning,list_snps[i].number_in_line)
                    temp_snp.file = list_snps[i].file
                    temp_snp.isolate = list_snps[i].isolate
                list_snps.append(temp_snp)
                
                ###WTF does this line even do? 
                ###seems like it just adds bases to the N possibilities and increases the numbers which results in increased occurrence
                ###uncertain bases should really just be removed
                list_replacements[-1].append(list_snps[i].wt)
            #for all other baes, 
            else:
                list_snps[i].alt = list_replacements[index_value][0] #replaces the ambiguous alternative with its first possibility
                #this loop adds the remaining alternatives to the list of SNPs
                for replacement in list_replacements[index_value][1:]:
                    temp_snp = SNP(list_snps[i].wt, replacement,list_snps[i].pos,'bogus/bogus',list_snps[i].line_beginning,list_snps[i].number_in_line)
                    temp_snp.file = list_snps[i].file
                    temp_snp.isolate = list_snps[i].isolate
                list_snps.append(temp_snp)
        i += 1