#This is the blast output reader
#note that the format must be outfmt 1 AKA Query-anchored showing identities

import os
import subprocess
import meta_pk as mpk
import time
import multiprocessing

#blast_filename = args.blast_filename
home_directory_path = '/storage/work/p/pcf5065/'
#this class is used to store the information from the ScanFilter results regarding pks
class putative_pk:
    def __init__(self, five_prime_boundary, three_prime_boundary, five_bond, three_bond, avg_z, median_z, high_z, num_pks, label, label_location, admin, length, nearby_labels, sequence):
        self.five_prime_boundary = five_prime_boundary
        self.actual_five_prime_boundary = five_prime_boundary
        self.three_prime_boundary = three_prime_boundary
        self.five_bond = five_bond
        self.three_bond = three_bond
        self.avg_z = avg_z
        self.median_z = median_z
        self.high_z = high_z
        self.num_pks = num_pks
        self.label = label
        self.label_location = label_location
        self.admin = admin
        self.length = length
        self.nearby_labels = nearby_labels
        self.sequence = sequence
        self.is_snp = False
        self.snp_string = ''
        self.relevant_sequence = sequence

    def create_SNP(self, SNP):
        self.sequence = rna_to_dna(self.sequence)
        print(self.sequence)
        print(self.five_prime_boundary, self.three_prime_boundary, self.length)
        print(int(SNP.pos)-int(self.five_prime_boundary))
        if self.sequence[int(SNP.pos)-int(self.five_prime_boundary)] == SNP.wt:
            snp_sequence = self.sequence[:SNP.pos-int(self.five_prime_boundary)] + '|' +SNP.alt +'|'+ self.sequence[SNP.pos-int(self.five_prime_boundary)+1:]
            #if SNP[2] < 61:
            #    snp_sequence = snp_sequence[:SNP[2]+60]
            #elif SNP[2] > len(snp_sequence)-1:
            #    snp_sequence = snp_sequence[SNP[2]-61:]
            #    SNP[2] = 61
            #else:
            #    snp_sequence = snp_sequence[SNP[2]-61:SNP[2]+60]
            #    SNP[2] = 61
            print(snp_sequence)
        else:
            print(SNP.extol(), 'WT base is not found in the sequence')
            return False
        snp_sequence = snp_sequence.replace('|', '')
        snp_pk = putative_pk(self.five_prime_boundary, self.three_prime_boundary, self.five_bond, self.three_bond, self.avg_z, self.median_z, self.high_z, self.num_pks, self.label, self.label_location, self.admin, self.length, self.nearby_labels, snp_sequence)
        snp_pk.is_snp = True
        return snp_pk

class SNP:
    def __init__(self, wt, alt, pos, isolate, line_beginning, number_in_line):
        self.wt = wt
        self.alt = alt
        self.pos = pos
        self.file = str(isolate)+'_blastn.txt'
        self.isolate = [isolate.split('/')[-1]]
        self.line_beginning = line_beginning
        self.number_in_line = number_in_line
        self.snp_dg = ''
        self.wt_dg = ''
        self.ddg = ''
        self.pk_results = 0
        self.rnasnp_results = 0
        self.occurence = 0
        self.snp_string = ''
        self.pk = 0
        self.ref_pk_filename = ''
        self.ref_snp_file_name = ''
    
    def initialize_list(self):
        self.list = [self.pos, self.wt, self.alt, self.occurence, self.snp_dg, self.wt_dg, self.ddg, self.pk_results, self.rnasnp_results]

    def extol(self):
        self.initialize_list()
        new_list = []
        for i in self.list:
            new_list.append(str(i))
        extol = '\t'.join(new_list)
        return extol

def get_snp_location_from_list(SNP):
    return int(SNP[1])
#this put blastn and snp outputs into foldesr
def organize_genome_folder():
    try:
        os.mkdir('blastn_outputs_Pass')
    except:
        pass
    try:
        os.mkdir('collected_snps_Pass')
    except:
        pass

    for file in os.listdir():
        if 'Pass' not in file:
            if 'blastn' in file:
                os.rename(file, 'blastn_outputs_Pass/'+file)
            if 'snp' in file:
                os.rename(file, 'collected_snps_Pass'+file)

def organize_genome_folder_SARS():
    try:
        os.mkdir('blastn_outputs_Pass')
    except:
        pass
    try:
        os.mkdir(home_directory_path + 'SARS_CoV_2/collected_snps_Pass')
    except:
        pass

    for file in os.listdir():
        if 'Pass' not in file:
            if 'blastn' in file:
                os.rename(file, 'blastn_outputs_Pass/'+file)
            if 'snp' in file:
                os.rename(file, home_directory_path + 'SARS_CoV_2/collected_snps_Pass/'+file)

#this function takes fasta files that have more than one genome in them and separates them into individual files
def extract_fasta_files(input_file):
    input_file_open = open(input_file, 'r')
    input_lines = input_file_open.readlines()
    name = ' '
    sequence = ''
    for line in input_lines:
        if '>' in line and name == ' ':
            name = line
        elif '>' not in line and name != ' ':
            sequence += line
        elif '>' in line and name != ' ':
            temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', 'w')
            temp_file.write(name+'\n'+sequence)
            temp_file.close()
            sequence = ''
            name = line
    os.remove(input_file)
    temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', 'w')
    temp_file.write(name+'\n'+sequence)
    temp_file.close()

def extract_fasta_files_species(input_file):
    input_file_open = open(input_file, 'r')
    input_lines = input_file_open.readlines()
    name = ' '
    sequence = ''
    for line in input_lines:
        if '>' in line and name == ' ':
            name = line
        elif '>' not in line and name != ' ':
            sequence += line
        elif '>' in line and name != ' ':
            temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', 'w')
            temp_file.write(name+'\n'+sequence)
            temp_file.close()
            try:
                print(name[1:4]+'_Pass')
                os.mkdir(name[1:4]+'_Pass')
            except:
                pass
            os.rename(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', name[1:4]+'_Pass'+'/'+name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta')
            sequence = ''
            name = line
    ###os.remove(input_file)
    temp_file = open(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', 'w')
    temp_file.write(name+'\n'+sequence)
    temp_file.close()
    try:
        print(name[1:4]+'_Pass')
        os.mkdir(name[1:4]+'_Pass')
    except:
        pass
    os.rename(name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta', name[1:4]+'_Pass'+'/'+name.replace('/','').replace('>','').replace('|','').replace('\n','').replace('-','').replace(' ','')[:51]+'.fasta')

def create_SNP(sequence, SNP):
    sequence = rna_to_dna(sequence)
    if sequence[61] == SNP[1]:
        snp_sequence = sequence[:61] + '|' +SNP[3] +'|'+ sequence[62:]
        print(snp_sequence)
    else:
        print(SNP, 'wt base is not found in the sequence')
        return False
    snp_sequence = snp_sequence.replace('|', '')
    return snp_sequence

def dna_to_rna(sequence):
    for number, char in enumerate(sequence):
        if char == 'T':
            sequence = sequence[:number]+'U'+sequence[number+1:]
    return sequence

def rna_to_dna(sequence):
    for number, char in enumerate(sequence):
        if char == 'U':
            sequence = sequence[:number]+'T'+sequence[number+1:]
    return sequence

def individual_blast(query, subject, output):
    subprocess.run('blastn -query %s -subject %s -outfmt 1 -max_hsps 1 -out %s' % (query, subject, output),shell = True)

#blastn -query -subject -outfmt 1 -out
#this function is the master function that blasts all genomes in comparison to the reference genome and runs all subsequent programs

def sequence_remaining_snps(snp):
    ref_file = home_directory_path+ 'SARS_CoV_2/NC_045512.2_SARS_CoV_2_ref_seq.fasta'
    dir_name = 'putative_SNitches/'
    ref = open(ref_file, 'r')
    ref_lines = ref.readlines()
    wt_sequence = ''
    for line in ref_lines[1:]:
        wt_sequence += line.replace('\n','')
    snp_string = snp.wt+str(snp.pos)+snp.alt
    temp_file_name =snp_string+'.fasta'
    temp_file = open(temp_file_name, 'w')
    if snp.alt != '-':
        temp_file.write('>rumor2'+snp_string+'\n'+wt_sequence[int(snp.pos)-60:int(snp.pos)] + snp.alt+wt_sequence[int(snp.pos)+1:int(snp.pos)+61])
    elif snp.alt == '-':
        print('224 deletion')
        print(snp_string)
        temp_file.write('>rumor2'+snp_string+'\n'+wt_sequence[int(snp.pos)-60:int(snp.pos)] +wt_sequence[int(snp.pos)+1:int(snp.pos)+61])
    temp_file.close()
    return temp_file_name

def sequence_remaining_wt_ref(snp):
    ref_file = home_directory_path+ 'SARS_CoV_2/NC_045512.2_SARS_CoV_2_ref_seq.fasta'
    dir_name = 'putative_SNitches/'
    ref = open(ref_file, 'r')
    ref_lines = ref.readlines()
    wt_sequence = ''
    for line in ref_lines[1:]:
        wt_sequence += line.replace('\n','')
    snp_string = snp.wt+str(snp.pos)+snp.alt
    temp_file_name = snp_string+'.fasta'
    temp_file = open(temp_file_name, 'w')
    temp_file.write('>rumor'+snp_string+'\n'+wt_sequence[int(snp.pos)-60:int(snp.pos)] + snp.wt+wt_sequence[int(snp.pos)+1:int(snp.pos)+61])
    temp_file.close()
    return temp_file_name

def run_blastn(putative_pks_input_file):
    #this is the principle functional of the RSC module
    print(time.ctime(time.time()), 'time')
    #this tracks the number of runs that have occurred for naming purposes
    number_runs = open('run_ID.txt', 'r')
    number_runs_lines = number_runs.readlines()
    print(number_runs_lines)
    print(number_runs_lines[0].replace('\n', ''), 'GGGG')
    num_prev_run = int(number_runs_lines[0].replace('\n', ''))
    this_run_num = num_prev_run + 1
    number_runs.close()
    number_runs = open('run_ID.txt', 'w')
    number_runs.write(str(this_run_num)+'\n')
    number_runs.close()
    
    #this is the working directory of the pipeline
    os.chdir('SARS_CoV_2')
    #The ref seq of SARS-CoV-2
    ref_file = '../NC_045512.2_SARS_CoV_2_ref_seq.fasta'
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
        #print(snp.extol(), 'final_check')
        #print(snp.isolate)
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


#this function takes ambiguous bases and formats them as separate possible SNPs in the list of SNPs
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

def snp_already_analyzed(strain_id, wt,loci,mu):
    strain = strain_id.split('/')[1].split('_')[0]

    #strain = strain_id.split('_')[0]
    #we're currently in the SARS_CoV_2 folder
    print(os.getcwd(),'407')
    os.chdir('collected_snps_Pass')
    for file in os.listdir(os.getcwd()):
        if strain in file:
            print('434',strain,strain_id,file)

            open_file = open(file,'r')
            lines = open_file.readlines()
            for line in lines[1:]:
                data = line.split(' ')
                if (wt == data[0]) and (mu == data[2]) and (int(loci) == int(data[1])):
                    print('already run, skipping:',str(line), strain)
                    os.chdir('../')
                    return True
    os.chdir('../')
    return False

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
                                            #list_location_snp.append([location, 'DEL '+item3, line_beginning_loci+location, 'del',line_beginning_loci])

        #this functions to determine what the mutant base was for
        #  the SNP
#        if addition == True:
 #           prev_items = file_lines[number_line-1].split(' ')
#            for item in prev_items:
#                if ('A' in item) or ('T' in item) or ('C' in item) or ('G' in item):
 #                   #for snpoly in temp_location_snp:
#                    snpoly = temp_location_snp
#                    for number, character in enumerate(item):
#                        if number == snpoly[0] and character != '.':
#                            list_location_snp.append([snpoly[0], snpoly[1], snpoly[2], character, snpoly[3]])
#                    break
    add_snps_nomenclature(list_location_snp)
    return list_location_snp

#this function formats the data calculated by ScanFilter
def extract_putative_pks(putative_pks_file_input):
    list_pks = []
    putative_pks_file = open(putative_pks_file_input, 'r')
    putative_pk_file_lines = putative_pks_file.readlines()
    for line in putative_pk_file_lines[1:]:
        items = line.split('\t')
        temp_pks = putative_pk(items[0],items[1],items[2],items[3],items[4],items[5],items[6],items[7],items[8],items[9],items[10],items[11],items[12],items[13].replace('\n',''))
        list_pks.append(temp_pks)
    return list_pks

#this makes a sequence file in fasta format
def make_struc_fasta(ID, snp_string, sequence):
    snp_name = snp_string
    temp_fasta_name = str(ID)+snp_name+'.txt'
    temp_fasta = open(temp_fasta_name, 'w')
    temp_fasta.write('>'+str(ID)+snp_name+'\n'+sequence)
    temp_fasta.close()
    temp_fasta = open(temp_fasta_name, 'r')
    return temp_fasta_name

#OBSOLETE function, counts number of pseudoknot base pairs
def nested_pairs(file_lines):
    pseudoknot_count = 0
    smallest_num = 0
    largest_num = 0
    running_small = 0
    running_large = 0
    adjacent = False
    known_bases = []
    for line in file_lines[1:]:
        indiv_elements = line.split(' ')
        indiv_elements[2].replace('\n','')
        indiv_elements[0] = int(indiv_elements[0])
        indiv_elements[2] = int(indiv_elements[2])
        if (indiv_elements[0] not in known_bases) and (indiv_elements[2] != 0):
            known_bases.append(indiv_elements[0])
            known_bases.append(indiv_elements[2])
            if (smallest_num == largest_num == running_small == running_small) or (indiv_elements[0]>largest_num):
                smallest_num = running_small = indiv_elements[0]
                largest_num = running_large = indiv_elements[2]
            elif (indiv_elements[0] in range(running_small, running_large+1)) and (indiv_elements[2] in range(running_small, running_large+1)):
                running_small = indiv_elements[0]
                running_large = indiv_elements[2]
            elif ((indiv_elements[0] in range(smallest_num, largest_num+1)) and (indiv_elements[2] in range(smallest_num, largest_num+1)) and (indiv_elements[0] not in range(running_small, running_large+1)) and (indiv_elements[2] not in range(running_small, running_large+1))):
                running_small = indiv_elements[0]
                running_large = indiv_elements[2]
            elif (adjacent == False):
                pseudoknot_count += 1
                adjacent = True      
                continue  
        adjacent = False
    return pseudoknot_count

def run_spot_wt_and_snp(ID, ID2, snp_string, snp_sequence, wt_sequence, input_file, working_dir = ''):
    wt_ID = 'ref_wt_'+ID
    try:
        os.mkdir('putative_SNitches')
    except:
        pass
    os.chdir('putative_SNitches')
    try:
        os.mkdir('SPOT-RNA_pseudoSNitch_data')
    except:
        pass
    os.chdir('SPOT-RNA_pseudoSNitch_data')
    try:
        os.mkdir(str(ID)+'.drct')
    except:
        pass
    os.chdir(str(ID)+'.drct')
    temp_snp_fasta = make_struc_fasta(ID, snp_string, snp_sequence)
    temp_wt_fasta = make_struc_fasta(wt_ID, '', wt_sequence)
    print(os.getcwd(), working_dir)
    os.chdir('../../../')
    print(os.getcwd())
    if working_dir == '':
        os.chdir('programs/SPOT-RNA')
    else:
        os.chdir('../programs/SPOT-RNA')
    spot_rna_output_path= '../../'+working_dir+'putative_SNitches/SPOT-RNA_pseudoSNitch_data/'+str(ID)+'.drct/'
    spot_rna_output_path_huh= '../../'+working_dir+'putative_SNitches/SPOT-RNA_pseudoSNitch_data/'+str(ID)+'.drct/'
    subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+str(temp_wt_fasta)+' --outputs '+spot_rna_output_path_huh, shell=True)
    subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+str(temp_snp_fasta)+' --outputs '+spot_rna_output_path_huh, shell=True)
    os.chdir(spot_rna_output_path)
    comparison = mpk.scorer_for_pks(ID, str(ID)+snp_string+'.ct', 'ref_wt_'+ID+'.ct', 2,2)
    os.chdir('../../')
    track = store_real_pseudoSNitches(ID, snp_string, comparison, 'spot_rna', input_file)
    return [comparison, track]

#this function runs SPOT-RNA
def run_spot_snp(ID, ID2, snp_string, sequence, degree_support, input_file):
    #based on the function run_spot in ScanFold-Umbrella
    #Changed to suit the purposes of this script
    #PATH
    try:
        os.mkdir('SPOT-RNA_pseudoSNitch_data')
    except:
        pass
    os.chdir('SPOT-RNA_pseudoSNitch_data')
    try:
        os.mkdir(str(ID)+'.drct')
    except:
        pass
    os.chdir(str(ID)+'.drct')
    temp_fasta = make_struc_fasta(ID, snp_string, sequence)
    #PATH
    print(os.getcwd(), 'DDD')
    os.chdir('../../../../')
    print(os.getcwd(), 'EE')
    os.chdir('programs/SPOT-RNA')
    spot_rna_output_path= '../../SARS_CoV_2/putative_SNitches/SPOT-RNA_pseudoSNitch_data/'+str(ID)+'.drct/'
    spot_rna_output_path_huh= spot_rna_output_path
    subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+str(temp_fasta)+' --outputs '+spot_rna_output_path_huh, shell=True)
    os.chdir(spot_rna_output_path)
    #num_pks = find_num_pks(str(pk.five_prime_boundary)+pk.snp_string+'.bpseq')
    #this runs the pseudoknot checker and analyzes the prediction
    is_spot_sig = SPOT_RNA_comparison_pk(ID, ID2, snp_string)
    #pk.num_pks = num_pks
    #PATH
    os.chdir('../../')
    results = store_real_pseudoSNitches(ID, snp_string, is_spot_sig, 'spot_rna', input_file)
    return results

#this function serves to compare the predicted structure of the SNP with the known from ScanFilter
def SPOT_RNA_comparison_pk(ID, ID2, snp_string):
    comparison = mpk.scorer_for_pks(ID, str(ID)+snp_string+'.ct', home_directory_path+'known_spot_data/'+str(ID2)+".drct/"+str(ID2)+'.ct', 2,2)
    return comparison

#this function sets aside a sequence when it is found that there are mulptiple SNPs in a region of interest
def multiple_snps_in_pk_region(input_file, SNPs):
    print(input_file)
    real_name = input_file.split('/')
    print(real_name)
    print(os.getcwd())
    #PATH
    os.chdir(real_name[0])
    try:
        os.mkdir('Pass_multiple_snps_in_pk_region')
    except:
        pass
    os.chdir('Pass_multiple_snps_in_pk_region')
    multiple_snps_file = open('mulitiple_snp_strains.txt', 'a')
    output_str = ''
    for snp in SNPs:
        output_str += snp.extol()
    multiple_snps_file.write(real_name[1]+' '+str(output_str)+'\n')
    os.chdir('../../')

#this function sets aside regions that have an INDEL and an SNP in them
def indel(input_file, SNP):
    #PATH
    real_name = input_file.split('/')
    os.chdir(real_name[0]+'/')
    try:
        os.mkdir('Pass_sequences_with_indel_in_putative_pk_region')
    except:
        pass
    os.chdir('Pass_sequences_with_indel_in_putative_pk_region')
    indel_file = open('indels_in_regions_of_interest','w')
    indel_file.write(real_name[1] + ' '+ str(SNP)+'\n')
    #os.system('cp %s Pass_sequences_with_indel_in_putative_pk_region/%s' % (real_name[1], real_name[1]))

    #os.rename(real_name[1], 'Pass_sequences_with_indel_in_putative_pk_region/'+real_name[1])
    #non_blast = real_name[1].replace('_blastn', '')
    #os.system('cp %s Pass_sequences_with_indel_in_putative_pk_region/%s' % (real_name[1], real_name[1]))

    #os.rename(non_blast, 'Pass_sequences_with_indel_in_putative_pk_region/'+non_blast)
    os.chdir('../../')

def run_RNAstructure_fold(sequence_fasta, snp):
    fold_output_name = 'Fold_output.ct'
    fold_output = open(fold_output_name, 'w')
    subprocess.run('Fold -mfe %s %s' % (sequence_fasta, fold_output_name),shell=True)
    efn2_output_name = 'efn2_output.txt'
    efn2_output = open(efn2_output_name, 'w')
    subprocess.run('efn2 %s %s' % (fold_output_name, efn2_output_name), shell=True)
    efn2_reader = open(efn2_output_name, 'r', encoding='utf-8')
    mfe_lines = efn2_reader.readlines()
    info_string = mfe_lines[0].split('= ')[1].replace('\n','').replace('\xb1', '+/-')
    os.remove(fold_output_name)
    os.remove(efn2_output_name)
    return info_string


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

def remove_excess():
    for file in os.listdir(os.getcwd()):
        if '_SNP' in file or 'RNAsnp' in file:
            os.remove(file)

#this makes the sequence files that the following programs use
def make_seq_files(pk):
    #these files are not in fasta format
    if '.fasta' not in str(pk.five_prime_boundary):
        temp_fasta_name = str(pk.five_prime_boundary)+'.fasta'
    else:
        temp_fasta_name = str(pk.five_prime_boundary)
    if pk.is_snp == True:
        temp_fasta_name = temp_fasta_name.replace('.fasta', pk.snp_string+'_SNP.fasta')
    temp_fasta = open(temp_fasta_name, 'w')
    pk.relevant_sequence = dna_to_rna(pk.relevant_sequence)
    #temp_fasta.write('>'+str(temp_fasta_name).replace('.fasta', '')+'\n')
    temp_fasta.write(pk.relevant_sequence)
    temp_fasta.close()
    return temp_fasta_name    

#reads the data from the RNAsnp output and determines if there have been a change
def RNAsnp_significance(RNAsnp_output, RNAsnp_output_name, degree_support):
    d_max = False
    r_min = False
    change = False
    RNAsnp_output.close()
    RNAsnp_output = open(RNAsnp_output_name, 'r')
    RNAsnp_output_lines = RNAsnp_output.readlines()
    print(RNAsnp_output_lines)
    for num_line, line in enumerate(RNAsnp_output_lines):
        split_lines = line.split('\t')
        if 'SNP' in split_lines[0]:
            try:
                data_line = RNAsnp_output_lines[num_line+1].split('\t')
                break
            except:
                error_RNAsnp_file = open('RNAsnp_error_log.txt', 'a')
                error_RNAsnp_file.write(RNAsnp_output_name+'\n')
                return None
    if len(data_line) == 1:
        error_RNAsnp_file = open('RNAsnp_error_log.txt', 'a')
        error_RNAsnp_file.write(RNAsnp_output_name+'\n')
        return None
    if float(data_line[6]) <= 0.1:
        d_max = True
        degree_support+=1
    if float(data_line[9].replace('\n', '')) <= 0.1:
        r_min = True
        degree_support+=1 
    if r_min == True or d_max == True:
        change = True
    output_string = ''
    for line in RNAsnp_output_lines:
        output_string += line
    return [change, degree_support, 'RNAsnp', output_string]

#reads the data from the SNPfold output and reports it
def SNPfold_significance(SNPfold_output_name, degree_support):
    cc_bpprob = 0.0
    cc_shannon = 0.0
    bp_p = 0.0
    shannon_p = 0.0
    change = False
    degree_support = 0
    #SNPfold_output.close()
    SNPfold_output = open(SNPfold_output_name, 'r')
    SNPfold_output_lines = SNPfold_output.readlines()
    data_line = SNPfold_output_lines[1].split('\t')
    print(data_line)
    cc_bpprob = float(data_line[1])
    bp_p = float(data_line[2])
    cc_shannon = float(data_line[3])
    shannon_p = float(data_line[4])
    if cc_bpprob <= 0.8:
        change = True
        degree_support += 1
    if bp_p <= 0.1 and change == True:
        degree_support +=1
    if cc_shannon <= 0.8:
        change = True
        degree_support += 1
        if shannon_p < 0.1:
            degree_support += 1
    print (cc_bpprob, bp_p,cc_shannon, shannon_p)
    return [change, degree_support, 'SNPfold' ]
    print (SNPfold_output_lines)

#this reads the data form SNPfold output if the -A option is off
def SNPfold_significance_multiple_SNPs(SNPfold_output_name, degree_support):
    cc_bpprob = 0.0
    cc_shannon = 0.0
    change = False
    degree_support = 0
    #SNPfold_output.close()
    SNPfold_output = open(SNPfold_output_name, 'r')
    SNPfold_output_lines = SNPfold_output.readlines()
    print(SNPfold_output_lines, 'JJJJJ')
    if 'ERROR' in SNPfold_output_lines[0]:
        return [change, degree_support, 'SNPfold' ]
    data_line = SNPfold_output_lines[1].split('\t')
    cc_bpprob = float(data_line[1])
    cc_shannon = float(data_line[2])
    if cc_bpprob <= 0.8:
        change = True
        degree_support += 1
    if cc_shannon <= 0.8:
        change = True
        degree_support += 1
    print (cc_bpprob, cc_shannon)
    return [change, degree_support, 'SNPfold' ]
    print (SNPfold_output_lines)

def remuRNA_significance(remuRNA_output_name):
    change = False
    output_string = ''
    degree_support = 0
    remuRNA_output = open(remuRNA_output_name, 'r')
    remu_lines = remuRNA_output.readlines()
    for line in remu_lines:
        output_string += line
    tab_items = remu_lines[1].split('\t')
    H = float(tab_items[4])
    if H < .8:
        change = True
        degree_support = 1
    return [change, degree_support,'remuRNA', output_string]

def run_RNAsnp(name, snp_string, wt_filename, snp_SNP_filename, degree_support, input_file, dir_name, default=True):
    if default == True:
        input_ID = input_file.split('/')[1][:-11]
    else:
        input_ID = input_file
    RNAsnp_output_name = 'RNAsnp_output%s_%s.txt' % (snp_string, input_ID)
    RNAsnp_output = open(RNAsnp_output_name, 'w')
    print(wt_filename, snp_SNP_filename, 'VVVVV')
    subprocess.run('RNAsnp -f %s -s %s -w 100 > %s' % (wt_filename, snp_SNP_filename, RNAsnp_output_name), shell=True)
    is_RNAsnp_sig = RNAsnp_significance(RNAsnp_output, RNAsnp_output_name, degree_support)
    if is_RNAsnp_sig != None:
        print(is_RNAsnp_sig[:3], '\n')
        store_real_riboSNitches(name, snp_string, is_RNAsnp_sig, is_RNAsnp_sig[3], input_file, default)
        os.rename(RNAsnp_output_name,dir_name+'/'+RNAsnp_output_name )
    return is_RNAsnp_sig

def run_remuRNA(snp_pk_five_prime_boundary, snp_string, pk_filename, snp_SNP_filename, degree_support, input_file, dir_name):
    input_ID = input_file.split('/')[1][:-11]
    remuRNA_output_name = 'remuRNA_output%s_%s.txt' % (snp_string, input_ID)
    remuRNA_output = open(remuRNA_output_name, 'w')
    pk_output = open(pk_filename, 'r')
    pk_seq = pk_output.readlines()
    remu_fasta_format = open('remuRNA_input.fasta', 'w')
    remu_fasta_format.write('>%s\n%s*%s' % (input_ID, pk_seq, snp_string))
    remu_fasta_format.close()
    os.chdir('../programs/remuRNA/data')
    subprocess.run('remuRNA %s > %s' % ('../../../putative_SNitches/'+remu_fasta_format, '../../../putative_SNitches/'+remuRNA_output_name), shell=True)
    os.chdir('../../../putative_SNitches')
    os.remove('remuRNA_input.fasta')
    is_remuRNA_sig = remuRNA_significance(remuRNA_output_name)
    print(is_remuRNA_sig[:3], '\n')
    store_real_riboSNitches(snp_pk_five_prime_boundary, snp_string, is_remuRNA_sig, is_remuRNA_sig[3], input_file)
    os.rename(remuRNA_output_name, dir_name+'/'+remuRNA_output_name)
    return is_remuRNA_sig

#this function sets up the intermediate file for ACI so that SNPfold can run
def set_up_intermediate_SNPfold(ID, snp_string, wt_filename, input_file, degree_support,this_run_num, snp_list, dir_name, working_dir = '', default=True):
    intermediate_snpfold_storage_file = open('intermediate.txt', 'a')
    intermediate_snpfold_storage_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n '%(str((ID)), str(snp_string), wt_filename, input_file, str(0), str(this_run_num), str(snp_list), dir_name, str(default)))
    intermediate_snpfold_storage_file.close()

#this runs both of the riboSNitch-finding programs to determine if there is a change in the secondary structure due to a SNP, or multiple SNPs
def riboSNitch_programs(pk, snp_pk, SNP, degree_support, input_file, this_run_num):
    snp_string = snp_pk.snp_string.replace('_', '-')
    pk_filename = make_seq_files(pk)
    snp_pk_filename = make_seq_files(snp_pk)
    print(snp_pk_filename, 'KKKKK', os.getcwd())
    snp_SNP_filename = snp_string+'.fasta'
    print(snp_string)
    snp_SNP_file = open(snp_SNP_filename, 'w')
    snp_SNP_file.write(snp_string)
    snp_SNP_file.close()
    dir_name = 'riboSNitch_programs_data'
    try:
        os.mkdir(dir_name)
    except:
        pass
    print('Analysing with RNAsnp...')
    is_rnasnp_sig = run_RNAsnp(snp_pk.five_prime_boundary, snp_string, pk_filename, snp_SNP_filename, degree_support, input_file, dir_name)
    
    os.chdir('../')
    ref_pk_filename = 'ref_'+str(pk.five_prime_boundary)+snp_string+'.fasta'
    ref_snp_file_name = open(ref_pk_filename, 'w')
    SNP.ref_pk_filename = ref_pk_filename
    SNP.ref_snp_file_name = ref_snp_file_name
    ref_snp_file_name.write('>'+ref_pk_filename+'\n'+pk.relevant_sequence)
    os.chdir('putative_SNitches')
#    intermediate_snpfold_storage_file = open('intermediate.txt', 'a')
#    try:
#        intermediate_snpfold_storage_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n '%(str(snp_pk.five_prime_boundary), str(snp_string), ref_pk_filename, input_file, str(is_rnasnp_sig[1]), this_run_num, str(SNP),dir_name, 'True'))
#    except:
#        intermediate_snpfold_storage_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n '%(str(snp_pk.five_prime_boundary), str(snp_string), ref_pk_filename, input_file, str(is_rnasnp_sig), this_run_num, str(SNP),dir_name, 'True'))
#      
#    intermediate_snpfold_storage_file.close()
    
    os.rename(snp_string+'.fasta', dir_name+'/'+snp_string+'.fasta')

    if is_rnasnp_sig != None:
        return is_rnasnp_sig[1]
    else:
        return None

def run_probknot(pk, snp_pk, SNP, input_file):
    try:
        os.mkdir('prob_knot_verification')
    except:
        pass
    os.chdir('prob_knot_verification')
    probknot_output_ct = str(pk.five_prime_boundary)+pk.snp_string+'_probknot.ct'
    probknot_snp_output_ct = str(pk.five_prime_boundary)+snp_pk.snp_string+'_snp_probknot.ct'
    scorer_output = str(pk.five_prime_boundary)+snp_pk.snp_string+'_scorer_output.txt'

    pk_seq_file = '../'+str(pk.five_prime_boundary)+'.fasta'
    print(pk_seq_file, 'FFFFF', os.getcwd())
    subprocess.run('ProbKnot --sequence %s %s' % (pk_seq_file, probknot_output_ct), shell=True)
    snp_seq_file = pk_seq_file.replace('.fasta', snp_pk.snp_string+'_SNP.fasta')
    
    ###this eliminates an error that arises when dels (-) are in the sequence for probknot
    edit_out_dels = open(snp_seq_file,'r')
    old_lines = edit_out_dels.readlines()
    old_line = old_lines[0].replace('-','')
    edit_out_dels.close()
    new_file = open(snp_seq_file,'w')
    new_file.write(old_line)
    new_file.close()

    print(snp_pk.snp_string, snp_seq_file, 'GGGGG')
    subprocess.run('ProbKnot --sequence %s %s' % (snp_seq_file, probknot_snp_output_ct), shell=True)

    is_prob_knot_sig = mpk.scorer_for_pks(pk.five_prime_boundary, probknot_snp_output_ct, probknot_output_ct)
    #PATH
    os.chdir('../')
    results = store_real_pseudoSNitches(snp_pk.five_prime_boundary, snp_pk.snp_string, is_prob_knot_sig, 'probknot', input_file)
    return results

#this is the more general version of probknot
def run_probknot_general(ID, snp_string, snp_filename, wt_filename, input_file):
    wt_ID = 'ref_wt_'+ID
    try:
        os.mkdir('prob_knot_verification')
    except:
        pass
    os.chdir('prob_knot_verification')
    probknot_output_ct = wt_ID+'_probknot.ct'
    probknot_snp_output_ct = ID+snp_string+'_snp_probknot.ct'
    scorer_output =ID+snp_string+'_scorer_output.txt'

    wt_fasta_name = '../../'+wt_filename
    snp_fasta_name = '../../'+snp_filename

    subprocess.run('ProbKnot --sequence %s %s' % (wt_fasta_name, probknot_output_ct), shell=True)
    subprocess.run('ProbKnot --sequence %s %s' % (snp_fasta_name, probknot_snp_output_ct), shell=True)

    is_prob_knot_sig = mpk.scorer_for_pks(ID, probknot_snp_output_ct, probknot_output_ct)
    os.chdir('../')
    results = store_real_pseudoSNitches(ID, snp_string, is_prob_knot_sig, 'probknot', input_file)
    return results

#stores any SNP that has been found to alter structure of the pseudoknot
def store_real_pseudoSNitches(ID, snp_string, comparison, test_name, input_file):
    if comparison[4] == comparison[5] and comparison[6] == comparison[7]:
        pass
        return None
    else:
        try:
            os.mkdir('putative_pseudoSNitches')
        except:
            pass
        #os.chdir('putative_pseudoSNitches')
        if comparison[5] == 0:
            output_string = 'Pseudoknot PPV: %s / %s = %s\n' % (str(comparison[6]), str(comparison[7]), str(float(comparison[6]/float(comparison[7]))))
        elif comparison[7] == 0:
            output_string = 'Pseudoknot Sensitivity: %s / %s = %s\n' % (str(comparison[4]), str(comparison[5]), str(float(comparison[4]/float(comparison[5]))))
        else:
            output_string = 'Pseudoknot Sensitivity: %s / %s = %s\nPseudoknot PPV: %s / %s = %s' % (str(comparison[4]), str(comparison[5]), str(float(comparison[4])/float(comparison[5])),str(comparison[6]), str(comparison[7]), str(float(comparison[6]/float(comparison[7]))))
        if int(comparison[5]) < int(comparison[4]):
            track = 'Strengthened'
        else:
            track = 'Weakened'
        print(ID)
        print(snp_string)
        print(test_name)
        print(track)
        print('putative_pseudoSNitches/'+ID+snp_string+'_'+test_name+'_'+track+'.txt')
        temp_file = open('putative_pseudoSNitches/'+ID+snp_string+'_'+test_name+'_'+track+'.txt', 'a')
        temp_file.write(ID+'_'+snp_string+'\n'+test_name+'\n'+input_file+'\n'+output_string+'\n'+track)
        print('The SNP causes a predicted PseudoSNitch\n')
        return track

#sotre any SNP that causes a riboSNtich          
def store_real_riboSNitches(name,pk_snp_string, validation, data, input_file, default=True):
    name = str(name).replace(' ', '')
    if default == True:
        input_ID = input_file.split('/')[1][:-11]
    else:
        input_ID = input_file
    if validation[0] == True:
        try:
            os.mkdir('putative_riboSNitches')
        except:
            pass
        temp_file = open('putative_riboSNitches/'+str(name)+pk_snp_string+'_'+validation[2]+'.txt', 'a')
        temp_file.write('>'+str(name)+'_'+pk_snp_string+'\n'+validation[2]+'\n'+input_file+'\n'+data)
        print('The SNP causes a predicted RiboSNitch\n')

def store_rando_files(all_files):
    #os.chdir('putative_SNitches')
    SNP_file_name = 'SNP_files_Pass.drct'
    actual_seq_file = 'seq_files_Pass.drct'
    try:
        os.mkdir(SNP_file_name)
    except:
        pass
    try:
        os.mkdir(actual_seq_file)
    except:
        pass
    for file in os.listdir():
        if '.drct' not in file:
            if 'SNP' in file:
                os.rename(file, SNP_file_name+'/'+file)
            #elif '.fasta' in file:
            #    os.rename(file, actual_seq_file+'/'+file)
    for file in all_files:
        os.rename(file, actual_seq_file+'/'+file)
#check_snps_in_pks(blast_filename, pk_filename)
#probknot_pk_verification(pk_filename)

#run_blastn('putative_pks_SARS_CoV_2_ref_seq')

# A 36 T Occurence dgSNP dgWT ddg rbSN psSN