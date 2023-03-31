from ast import excepthandler
import multiprocessing as mp
from re import search
import time
import os
import argparse
import meta_pk as mpk
import sys

time_id = (time.ctime(time.time())).replace(' ', '_').replace(':', '.')

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', type=str, help='input reference sequence')
    parser.add_argument('-s', type=str, default='',help='input RSC library search parameter to create a library. No input requires r option.')
    parser.add_argument('-d', type=str, default=1000,help='input desired library depth. default 1000. Not required for l option')
    parser.add_argument('-c', type=str, help='input candidate regions in standard format')
    parser.add_argument('-g', type=str, help='gff3 file of the species (optional)')
    parser.add_argument('-l', type=str, default='', help='optional: input path to collected sequence library. Default: RSC module creates its own sequence library from an NCBI search (requires s option).')

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
    os.chdir(home_dir)
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
    os.system('esearch -db nucleotide -query "%s complete" | efetch -stop %s -format fasta > library_%s.fas'%(search, str(library_depth), run_ID))
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
    os.system('blastn -query %s -subject %s -outfmt 1 -max_hsps 1 -out %s' % (query, subject, output))

def run_blastn_on_lib(ref_seq, path_to_RSC_library, lib_depth):
    file_list = []
    orig_dir = os.getcwd()
    os.chdir(path_to_RSC_library)
    for numfile, file in enumerate(os.listdir(os.getcwd())):
        if numfile > int(lib_depth):
            break
        if '.fasta' in file:
            file_list.append(file)
    pool_manager_blastn(file_list)
    os.chdir(orig_dir)

def run_blastn_on_sequence(file):
    blastn_output_filename = file.replace('.fasta', '_blastn.txt')
    blastn_output_file = open(blastn_output_filename,'w')
    individual_blast(file, home_dir+ref_seq, blastn_output_filename)

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
            #accomodates insertions and adjusts snp loci based on that
            for num_place, place in enumerate(file_lines[number_line-1]):
                if place in ['A', 'G', 'C', 'T']:
                    start_seq_place = num_place
                    break
            #this for loop goes thorguh each entry in the data line
            for num_item, item in enumerate(line_items):
                #this defines where in the genome the SNP is located
                if str(item).isdigit() == True:
                    line_beginning_loci = int(item)
                    continue
                #this finds and identifies the wt base
                if '.' in item:
                    adjustment_for_ins = 0
                    for wt_number, wt_character in enumerate(item):
                        if wt_character == '-':
                            adjustment_for_ins += 1
                        if wt_character.capitalize() in ['A', 'G', 'C', 'T']:
                            adjustment_for_dels = 0
                            if len(file_lines[number_line+1]) > 4:
                                for top_lvl_num, possible_ins in enumerate(file_lines[number_line+1][:start_seq_place+wt_number+1]):
                                    if top_lvl_num == start_seq_place:
                                        continue
                                    if possible_ins == '\\':
                                        search_down = 2
                                        while (file_lines[number_line+search_down][top_lvl_num] == '|'):
                                            if file_lines[number_line+search_down+1][top_lvl_num] in ['A', 'G', 'C', 'T']:
                                                adjustment_for_dels += 1
                                                search_right = 1
                                                while file_lines[number_line+search_down+1][top_lvl_num+search_right] in ['A', 'G', 'C', 'T']:
                                                    adjustment_for_dels += 1
                                                    search_right += 1
                                                    if (search_right+top_lvl_num) >= len(file_lines[number_line+search_down+1]):
                                                        break
                                                search_left= 1
                                                while file_lines[number_line+search_down+1][top_lvl_num-search_left] in ['A', 'G', 'C', 'T']:
                                                    adjustment_for_dels += 1
                                                    search_left +=1
                                            if len(file_lines[number_line+search_down+1]) < 4:
                                                break
                                            search_down += 1
                            snp_loci = wt_number + line_beginning_loci-adjustment_for_ins+ adjustment_for_dels
                            prev_items = file_lines[number_line-1].split(' ')
                            #print(filename)
                            #print(snp_loci, wt_number)
                            for pitem in prev_items:
                                if ('A' in pitem) or ('T' in pitem) or ('C' in pitem) or ('G' in pitem):
                                    #for alt_number, alt_character in enumerate(pitem):
                                        #if (alt_number == wt_number) and (alt_character != '.'):
                                    alt_character = pitem[wt_number]
                                    #print(filename)
                                    #print(wt_character, snp_loci, alt_character)
                                    #print(line_beginning_loci)
                                    #print(wt_number)
                                    #print(adjustment_for_dels)
                                    #print(adjustment_for_ins)
                                    #print('---\n')
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
    os.chdir(home_dir)
    return path_to_RSC_lib+output_filename

def count_occurrences(lines):
    indecies_2_remove = []
    new_lines = []
    for n, line in enumerate(lines):
        if n not in indecies_2_remove:
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
        elif n in indecies_2_remove:
            occ = 1
            info = line.split('\t')
            acc = info[3].replace('\n','')
        new_line = '\t'.join(info[:3])+'\t'+str(occ)+'\t'+acc
        new_lines.append(new_line)
    indecies_2_remove.sort()
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
                    snp.region = int(cand.start)
                    snp.region_stop = int(cand.stop)
                    if check_if_snp_in_list(snp, list_cand_snps) == False:
                        list_cand_snps.append(snp)
    return [list_cand_snps, list_candidates]

def check_if_snp_in_list(SNP, list_cand_snps):
    if len(list_cand_snps) == 0:
        return False
    for snp in list_cand_snps:
        if (snp.pos == SNP.pos) and (snp.wt == SNP.wt) and (snp.alt == SNP.alt):
            return True
    return False

def make_wt_and_alt_sequence(SNP):
    wt_filename = 'wt_'+SNP.ID+'.fasta'
    wt_file = open(wt_filename,'w')
    alt_filename = 'alt_'+SNP.ID+'.fasta'
    alt_file = open(alt_filename, 'w')
    if ref_seq_string[SNP.pos-1] == SNP.wt:
        five_prime_flanker = -61
        three_prime_flanker = 60 + int(SNP.pos)
        if SNP.pos <= 61:
            five_prime_flanker = SNP.pos-1
        if (int(SNP.pos)+60) > (len(ref_seq_string)-1):
            three_prime_flanker = (len(ref_seq_string)-1)
        ref_seq_window = ref_seq_string[int(SNP.pos+five_prime_flanker):three_prime_flanker]
        wt_file.write('>wt_'+SNP.ID+'\n'+ref_seq_window)
        alt_seq_window = ref_seq_string[int(SNP.pos)+five_prime_flanker:int(SNP.pos)-1] + SNP.alt+ref_seq_string[int(SNP.pos):three_prime_flanker]
        alt_file.write('>alt_'+SNP.ID+'\n'+alt_seq_window)
        SNP.wt_window = ref_seq_window
        SNP.alt_window = alt_seq_window
        print(wt_filename, SNP.wt_window)
        return [wt_filename, alt_filename]
    else:
        print('WT BASE DOES NOT MATCH SNP BASE\n')
        print(ref_seq_string[SNP.pos-3:SNP.pos+3])
        print(SNP.wt, SNP.alt, SNP.pos)

def get_individual_ddg(SNP):
    print(SNP.__dict__)
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
    os.system('Fold -q -mfe %s %s' % (seq_filename, fold_output_name))
    efn2_output_name = 'efn2_output_%s.txt' % (ID)
    efn2_output = open(efn2_output_name, 'w')
    os.system('efn2 -q %s %s' % (fold_output_name, efn2_output_name))
    efn2_reader = open(efn2_output_name, 'r', encoding='utf-8')
    #efn2_reader = open(efn2_output_name, 'r')
    mfe_lines = efn2_reader.readlines()
    try:
        info_string = mfe_lines[0].split('= ')[1].replace('\n','')
    except:
        print(seq_filename)
        sys.exit()
    os.remove(fold_output_name)
    os.remove(efn2_output_name)
    return info_string    

def get_ddgs(snps_in_candidate_regions, ref_seq_string):
    #establishes the number of available threads
    p = mp.Pool(mp.cpu_count())
    #p = mp.Pool(1)
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
    if SNP.pos < 61:
        snp_loc = str(SNP.pos)
    else:
        snp_loc = str(61)
    snp_file.write(SNP.wt+snp_loc+SNP.alt)
    snp_file.close()
    return snp_filename

def run_RNAsnp(SNP):
    filenames = make_wt_and_alt_sequence(SNP)
    wt_filename = filenames[0]
    alt_filename = filenames[1]

    RNAsnp_output_name = 'RNAsnp_output_%s.txt' % (SNP.ID)
    RNAsnp_output = open(RNAsnp_output_name, 'w')
    
    SNP_file = make_snp_file(SNP)
    os.system('RNAsnp -f %s -s %s -w 100 > %s' % (wt_filename, SNP_file, RNAsnp_output_name))
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
                os.remove(wt_filename)
                os.remove(alt_filename)
                os.remove(SNP_file)
                SNP.r_min = 'ERROR'
                SNP.d_max = 'ERROR'
                return SNP
    if len(data_line) == 1:
        error_RNAsnp_file = open('RNAsnp_error_log.txt', 'a')
        error_RNAsnp_file.write(RNAsnp_output_name+'\n')
        SNP.r_min = 'ERROR'
        SNP.d_max = 'ERROR'
        os.remove(RNAsnp_output_name)
        return SNP
    
    SNP.r_min = data_line[6].replace('\n','')
    SNP.d_max = data_line[9].replace('\n','')
    os.remove(RNAsnp_output_name)
    os.remove(wt_filename)
    os.remove(alt_filename)
    os.remove(SNP_file)
    return SNP

def make_region_seq_fasta(SNP):
    wt_filename = 'wt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'.fasta'
    wt_file = open(wt_filename,'w')
    alt_filename = 'alt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'.fasta'
    alt_file = open(alt_filename, 'w')
    if ref_seq_string[SNP.pos-1] == SNP.wt:
        ref_seq_window = ref_seq_string[SNP.region:SNP.region_stop+1]
        wt_file.write('>wt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'\n'+ref_seq_window)
        alt_seq_window = ref_seq_string[SNP.region:int(SNP.pos)-1] + SNP.alt+ref_seq_string[int(SNP.pos):SNP.region_stop+1]
        alt_file.write('>alt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'\n'+alt_seq_window)
        SNP.wt_window = ref_seq_window
        SNP.alt_window = alt_seq_window
        return [wt_filename, alt_filename]
    else:
        print('WT BASE DOES NOT MATCH SNP BASE')

def get_SPOT_RNA(snps_in_candidate_regions,path_to_SPOT_RNA):
    for SNP in snps_in_candidate_regions:
        orig_dir = os.getcwd()
        os.chdir(path_to_SPOT_RNA)
        filenames = make_region_seq_fasta(SNP)
        wt_filename = filenames[0]
        alt_filename = filenames[1]
        if 'wt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'.ct' not in os.listdir(os.getcwd()):
            os.system('python3 SPOT-RNA.py --inputs '+str(wt_filename)+' --outputs ./')
        os.system('python3 SPOT-RNA.py --inputs '+str(alt_filename)+' --outputs ./')
        comparison = mpk.scorer_for_pks(SNP.region, 'alt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'.ct', 'wt_region_'+str(SNP.region)+'_'+str(SNP.region_stop)+'.ct', 2,2)
        os.remove(wt_filename)
        os.remove(alt_filename)
        #os.remove(wt_filename.replace('.fasta', '.bpseq'))
        #os.remove(wt_filename.replace('.fasta', '.ct'))
        #os.remove(wt_filename.replace('.fasta', '.prob'))
        os.remove(alt_filename.replace('.fasta', '.bpseq'))
        os.remove(alt_filename.replace('.fasta', '.ct'))
        os.remove(alt_filename.replace('.fasta', '.prob'))
        os.chdir(orig_dir)
        SNP.pk_bps = comparison[5]
        try:
            SNP.pk_sensitivity = comparison[4]/comparison[5]
        except:
            SNP.pk_sensitivity = ''
        try:
            SNP.pk_ppv = comparison[6]/comparison[7]
        except:
            SNP.pk_ppv = ''
    return snps_in_candidate_regions

def set_up_intermediate_for_SNPfold(snps_in_candidate_regions,RSC_library_path, run_ID):
    snp_data_libraryname = RSC_library_path +'library_SNPs_data_'+run_ID+'.txt'
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
            outlist.append(str(info[col]).replace('\xb1','+/-'))
        outline = '\t'.join(outlist) +'\n'
        if outline!='\n':
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

def get_ProbKnot(snps_in_candidate_regions):
    #establishes the number of available threads
    #p = mp.Pool(mp.cpu_count())
    p=mp.Pool(1)
    #performs the ddG calculations in parallel
    snps_in_candidate_regions = p.map(run_probknot_general, snps_in_candidate_regions)
    #terminates the pool. This was essential when adding multiple computational batches per data set.
    p.terminate()
    return snps_in_candidate_regions

#this is the more general version of probknot
def run_probknot_general(SNP):
    filenames = make_wt_and_alt_sequence(SNP)
    wt_filename = filenames[0]
    alt_filename = filenames[1]
    probknot_wt_output_ct = SNP.ID+'_wt_probknot.ct'
    probknot_alt_output_ct = SNP.ID+'_alt_probknot.ct'

    os.system('ProbKnot --sequence %s %s' % (wt_filename, probknot_wt_output_ct))
    os.system('ProbKnot --sequence %s %s' % (alt_filename, probknot_alt_output_ct))

    comparison = mpk.scorer_for_pks(SNP.ID, probknot_alt_output_ct, probknot_wt_output_ct)
        
    SNP.pb_pk_bps = comparison[5]
    try:
        SNP.pb_pk_sensitivity = comparison[4]/comparison[5]
    except:
        SNP.pb_pk_sensitivity = ''
    try:
        SNP.pb_pk_ppv = comparison[6]/comparison[7]
    except:
        SNP.pb_pk_ppv = ''  
    os.remove(probknot_wt_output_ct)
    os.remove(probknot_alt_output_ct) 
    os.remove(wt_filename)
    os.remove(alt_filename)
    return SNP

def clean_spot_dir(list_snps, path_to_path_dir):
    orig_dir = os.getcwd()
    os.chdir(path_to_path_dir)
    for snp in list_snps:
        files = make_region_seq_fasta(snp)
        wt_filename = files[0]
        alt_filename = files[1]
        os.remove(alt_filename)
        os.remove(wt_filename)
        try:
            os.remove(wt_filename.replace('.fasta', '.bpseq'))
            os.remove(wt_filename.replace('.fasta', '.ct'))
            os.remove(wt_filename.replace('.fasta', '.prob'))
        except:
            pass
    os.chdir(orig_dir)

def master():
    args = get_args()
    search = args.s
    library_depth = args.d
    candidate_filename = args.c
    gff_filename = args.g
    global ref_seq
    ref_seq = args.r
    global ref_seq_string
    ref_seq_string = get_ref_seq_string(ref_seq)
    global home_dir
    home_dir = os.getcwd()+'/'
    if args.l != '':
        RSC_library_path = args.l
        run_ID = args.l.replace('/','')+'_'+time_id
        os.chdir(RSC_library_path)
        for file in os.listdir(os.getcwd()):
            if '.fas' in file:
                extract_fasta_files(file)
        os.chdir(home_dir)
    else:
        run_ID = search.replace(' ','_')+'_'+str(library_depth)+'_'+time_id
        RSC_library_path = construct_library(search,library_depth,run_ID)
    run_blastn_on_lib(ref_seq, RSC_library_path, int(library_depth))
    SNP_lib_path = collect_SNPs_to_lib(RSC_library_path, run_ID)
    lists_snps_cand = get_snps_in_candidate_regions(candidate_filename, SNP_lib_path)
    snps_in_candidate_regions = lists_snps_cand[0]
    cand_regions = lists_snps_cand[1]
    snps_in_candidate_regions = get_ddgs(snps_in_candidate_regions, ref_seq_string)
    snps_in_candidate_regions = get_ProbKnot(snps_in_candidate_regions)
    snps_in_candidate_regions = get_SPOT_RNA(snps_in_candidate_regions, '../programs/SPOT-RNA')
    clean_spot_dir(snps_in_candidate_regions,'../programs/SPOT-RNA')
    snps_in_candidate_regions = get_RNAsnp(snps_in_candidate_regions)
    global gff_obj_list
    gff_obj_list = get_gff_lines(gff_filename)
    snps_in_candidate_regions = label_AAs(snps_in_candidate_regions, gff_filename)
    intermediate_filename = set_up_intermediate_for_SNPfold(snps_in_candidate_regions,RSC_library_path, run_ID)
    return intermediate_filename

file = master()
print('###')
print(file)
sys.exit(file)