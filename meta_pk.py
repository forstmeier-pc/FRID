#psuedoknot comparison programs
#the purpose of this script is to test the various pseudoknot-finding programs

import subprocess
import argparse
import os
#from openpyxl import Workbook
#from openpyxl import load_workbook
#from blast_output_reader import *

class known_pk:
    def __init__(self):
        self.pkb_num = 'NOT EXTRACTED'
        self.definition = 'NOT EXTRACTED'
        self.organism = 'NOT EXTRACTED'
        self.type = 'NOT EXTRACTED'
        self.keywords = 'NOT EXTRACTED'
        self.accesion = 'NOT EXTRACTED'
        self.location = 'NOT EXTRACTED'
        self.supported = 'NOT EXTRACTED'
        self.sequence = 'NOT EXTRACTED'
        self.bracket = 'NOT EXTRACTED'

    def extol(self):
        string = self.pkb_num + '\t' + self.definition + '\t' + self.organism + '\t' + self.accesion +'\t'+self.location+ '\t' + self.type + '\t' + self.keywords + '\t' + self.supported + '\t' + self.sequence +'\t'+self.bracket
        return string

    def normalize_brackets(self):
        self.bracket = self.bracket.replace(':', '.')
        self.bracket = self.bracket.replace('[', '<')
        self.bracket = self.bracket.replace(']', '>')

#This function extracts data about known pseudoknots from PseudoBase
def extract_pseudobase(sub_html):
    admin_warning = False
    filename2 = 'tmp_file2.txt'
    subprocess.run('curl %s > %s' % (sub_html, filename2), shell=True)
    file2 = open(filename2, 'r')
    extraction_file_output_name = ''
    extracted_pseudoknot = known_pk()
    try:
        lines2 = file2.readlines()
    except:
        return
    #this loop graps important data from the curl extraction
    for line_num, line in enumerate(lines2):
        if 'PKB-number' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.pkb_num = lines2[line_num+2].replace('\n', '')
        if 'Definition' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.definition = lines2[line_num+2].replace('\n', '')
        if 'Organism' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.organism = lines2[line_num+2].replace('\n', '')
        if 'RNA&nbsp' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.type = lines2[line_num+2].replace('\n', '')
        if 'Keywords' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.keywords = lines2[line_num+2].replace('\n', '')
        if 'EMBL&' in line:
            #print (lines[line_num+2].split('"')[1][36:])
            extracted_pseudoknot.accesion = (lines2[line_num+2].split('"')[1][36:])
        if 'Supported' in line:
            #print (lines[line_num+2])
            extracted_pseudoknot.supported = lines2[line_num+2].replace('\n', '')
        if 'Bracket' in line:
            list_shifts = []
            for next_num, next_line in enumerate(lines2[line_num:]):
                if '%' in next_line:
                    list_shifts.append(next_num-1)
            shift = list_shifts[0]
            #print (lines[line_num+shift])
            #print (lines[line_num+shift+1])
            sequence = lines2[line_num+shift].replace('\n', '')
            sequence_index = sequence.split(' ')
            sequence_index = ' '.join(sequence_index).split()
            start_loci = sequence_index[1]
            try:
                end_loci = sequence_index[2].split('=')[1]
            except:
                try: 
                    end_loci = sequence_index[1].split('=')[1]
                except:
                    print ('\n'+'\n'+"ADMIN WARNING" +'\n'+'\n'+'NON-CONTINUOUS STRUCTURE, FILE NOT SAVED'+'\n'+'\n')
                    admin_warning = True
                    break
            extracted_pseudoknot.sequence=sequence_index[2].split('=')[0]
            bracket = lines2[line_num+shift+1].split(' ')
            bracket = ' '.join(bracket).split()[2].replace('\n', '')
            extracted_pseudoknot.bracket = bracket
            #this serves to interact with sequences and data that are on multiple lines
            for numb, jump in enumerate(list_shifts[1:]):
                shift = jump
                sequence_index2 = lines2[line_num+shift].replace('\n', '').split(' ')
                sequence_index2 = ' '.join(sequence_index2).split()
                bracket2 = lines2[line_num+shift+1].split(' ')
                try:
                    bracket2 = ' '.join(bracket2).split()[2].replace('\n', '')
                except:
                    bracket2 = ' '.join(bracket2).split()[1].replace('\n', '')
                #print(sequence_index2)
                if int(sequence_index2[1]) == int(end_loci)+1:
                    #print (lines[line_num+shift+1])  
                    end_loci = sequence_index2[2].split('=')[1]
                    #print(end_loci)
                    extracted_pseudoknot.sequence+=sequence_index2[2].split('=')[0]
                    extracted_pseudoknot.bracket += bracket2
                else:
                    print ('\n'+'\n'+"ADMIN WARNING" +'\n'+'\n'+'NON-CONTINUOUS STRUCTURE, FILE NOT SAVED'+'\n'+'\n')
                    admin_warning = True
                    break
            extracted_pseudoknot.location = start_loci+'-'+end_loci
            extracted_pseudoknot.normalize_brackets()

    #print(extracted_pseudoknot.extol())
    if admin_warning == False:
        store_pseudobase(extracted_pseudoknot)

#this function loops through all available links on the PKBASE page and extracts the data from them
def extract_all(master_html):
    http_start = 'http://www.ekevanbatenburg.nl/PKBASE/'
    filename = 'tmp_file.txt'
    subprocess.run('curl %s > %s' % (master_html, filename), shell=True)
    file = open(filename, 'r')
    extraction_file_output_name = ''
    lines = file.readlines()
    for line_num, line in enumerate(lines[90:]):
        if 'li' in line and '<a href' in line and '.HTML>' in line:
            print(line)
            line = line.split('=')
            html_indiv = http_start + line[1][:13]
            extract_pseudobase(html_indiv)

#this function stores the output from a pseudobase extraction
def store_pseudobase(pk):
    dir_name = 'pseudobase.drct'
    store_filename = 'stored_pseudobase_extractions_ALL.txt'
    os.chdir(dir_name)
    store_file = open(store_filename, 'a')
    q = pk.extol()
    store_file.write(q+'\n')
    store_file.close()
    os.chdir('../')

#This function uses the dot structure and the sequence extracted from pseudobase and builds a .dot file out of it
def make_dot_file(name, seq, dot_struc):
    temp_name = name+'.dot'
    temp_file = open(temp_name, 'a')
    temp_file.write('>'+name+'\n')
    temp_file.write(seq+'\n')
    temp_file.write(dot_struc)
    temp_file.close()
    return temp_name

#this function ensures that there are no duplicate files in the storage file
def check_for_doubles():
    os.chdir('pseudobase.drct')
    filename = 'stored_pseudobase_extractions_ALL.txt'
    file = open(filename, 'r')
    lines = file.readlines()
    known_names = []
    for line in lines:
        line = line.split('\t')
        name = line[0]
        if name not in known_names:
            known_names.append(name)
        else:
            print(name)

#this function converts the .dot files created from pseudobase information and converts
#it to a .ct file using the RNAstructure method
def pseudobase_2_ct():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('known_pk_dot_and_ct')
    except:
        pass
    filename = 'stored_pseudobase_extractions_ALL.txt'
    file = open(filename, 'r')
    lines = file.readlines()
    os.chdir('known_pk_dot_and_ct')
    for line in lines:
        line = line.split('\t')
        name = line[0]
        sequence = line[8]
        dot_struc = line[9]
        new_dot_file = make_dot_file(name, sequence, dot_struc)
        new_ct = name+'.ct'
        subprocess.run('dot2ct %s %s' % (new_dot_file, new_ct), shell=True)

#this function moves files from the base drct to other drcts
def clean_dir():
    os.chdir('pseudobase.drct')
    #os.mkdir('pairs')
    drct = os.listdir()
    print(os.getcwd())
    for num, file in enumerate(drct):
        file2 = file.split('.')
        if file2[0] == drct[num+1].split('.')[0]:
            os.rename(file, 'pairs/'+file)
            os.rename(drct[num+1], 'pairs/'+drct[num+1])

#this purpose of this function is to test if the pk_counter is consistent
def check_pk_counter():
    os.chdir('pseudobase.drct')
    os.chdir('known_pk_dot_and_ct')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        #fasta_file_name = make_fasta_files(name, sequence)
        ipknot_output_ct_name = name +'.ct'
        #subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, probknot_output_ct), shell=True)
        comparison = scorer_for_pks(name, name+'.ct', '../known_pk_dot_and_ct/')
        for place, value in enumerate(comparison):
            final_values[place] += value

#this function serves to compare the hits of ipknot to that of spot-rna to determine any overlap
#and to compare ipknot to the overlap of spot-rna with probknot
def compare_ip_to_spot():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('ipknot_massive')
    except:
        pass
    os.chdir('ipknot_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        #fasta_file_name = make_fasta_files(name, sequence)
        ipknot_output_ct_name = name +'.ct'
        #subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, probknot_output_ct), shell=True)
        comparison = scorer_for_pks(name, ipknot_output_ct_name, '../known_pk_dot_and_ct/')
        if comparison[6] == 0 and comparison[7] != 0:
            comparison = scorer_for_pks(name,'../probknot_massive/'+name+'.ct' , '',1)
            if comparison[6] == 0 and comparison[7] != 0:
                comparison = scorer_for_pks(name,'../spotrna_massive/'+name+'.drct/'+name+'.ct' , '',1,2)

                for place, value in enumerate(comparison):
                    final_values[place] += value
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values[0], final_values[1], str(float(final_values[0])/float(final_values[1]))))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values[2], final_values[3], str(float(final_values[2])/float(final_values[3]))))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values[4], final_values[5], str(float(final_values[4])/float(final_values[5]))))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values[6], final_values[7], str(float(final_values[6])/float(final_values[7]))))
    except:
        pass

#this function compared probknot to spot-rna to determine overlap
def compare_prob_to_spot(): 
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('probknot_massive')
    except:
        pass
    os.chdir('probknot_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    running_count = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        fasta_file_name = make_fasta_files(name, sequence)
        probknot_output_ct = name +'.ct'
        #subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, probknot_output_ct), shell=True)
        comparison = scorer_for_pks(name, probknot_output_ct, '../known_pk_dot_and_ct/')
        if comparison[6] == 0 and comparison[7] != 0:
            running_count += comparison[7]
            print(name, '@@@@')
            comparison = scorer_for_pks(name, '../spotrna_massive/'+name+'.drct/'+name+'.ct', '', 1,2)
        #number = export_to_excel('ProbKnot', name, comparison, number)
            for place, value in enumerate(comparison):
                final_values[place] += value
    print('RUNNING', running_count)
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values[0], final_values[1], str(float(final_values[0])/float(final_values[1]))))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values[2], final_values[3], str(float(final_values[2])/float(final_values[3]))))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values[4], final_values[5], str(float(final_values[4])/float(final_values[5]))))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values[6], final_values[7], str(float(final_values[6])/float(final_values[7]))))
    except:
        pass

def nested_checker_1bp_prediciton(accepted_pairs):
    nests = []
    pk_pairs = []
    known_bases = []
    end = False
    #this creates continuous 'nests' out of the accepted pairs
    for line, pair in enumerate(accepted_pairs):
        pair[0] = int(pair[0])
        pair[1] = int(pair[1])
        if pair[1] != 0 and pair[0] not in known_bases:
            if nests == []:
                temp_nest = nest_class()
                temp_nest.nest.append(pair)
                nests.append(temp_nest)
            elif pair[0] - 1 == accepted_pairs[line-1][0] and pair[1]+1 == accepted_pairs[line-1][1]:
                temp_nest.nest.append(pair)
            else:
                temp_nest = nest_class()
                temp_nest.nest.append(pair)
                nests.append(temp_nest)
            known_bases.append(pair[1])

    nests_to_remove = []

    #this removes nests that are totally self-nested, and therefore do not participate in pseudoknots
    for nest in nests:
        remove = True
        for pair in accepted_pairs:
            if pair[0] > nest.nest[-1][0]:
                #this conditional is the key... if the pair is outside of the nest and does not form a pk, then the nest is non-pseudoknotted
                if pair[1] != 0 and pair[0] <= nest.nest[-1][1]:
                    if is_nested(nest.nest[-1], pair) == False:
                        added_pair = [0,0]
                        if pair[1] < pair[0]:
                            added_pair[0] = pair[1]
                            added_pair[1] = pair[0]
                        else:
                            added_pair = pair
                        nest.competing_pairs.append(added_pair)
                        remove = False
        if remove == True:
            nests_to_remove.append(nest)
    for nest in nests_to_remove:
        nests.remove(nest)

    #this serves to combine nests that are separated by bulges or breaks and are nested
    for num_ses, suspect_nest in enumerate(nests):
        for nest_2 in nests[num_ses:]:
            same_comp = False
            if is_nested(suspect_nest.nest[-1], nest_2.nest[-1]) == True:
                for comp_base in suspect_nest.competing_pairs:
                    if comp_base in nest_2.competing_pairs:
                        same_comp = True
                if same_comp == True:
                    for pair_base in nest_2.nest:
                        suspect_nest.nest.append(pair_base)
                    nests.remove(nest_2)
    
    #this serves to combine two nests that are borken by the same pseudoknot
    if len(nests) % 2 != 0:
        for num_nest, nest in enumerate(nests):
            add_pairs = False
            for next_nest in nests[num_nest+1:]:
                for next_pair in next_nest.competing_pairs:
                    if next_pair in nest.competing_pairs: 
                        add_pairs = True
                        break
                if add_pairs == True:
                    for nest_next_pair in next_nest.nest:
                        nest.nest.append(nest_next_pair)
                    nests.remove(next_nest)

    #this is the final determiner of what combined nests of stuctures are pseudoknots and what are nested
    if len(nests) != 0 and len(nests) != 1:
        for num_nest, combined_nests in enumerate(nests):
            competing_struc = False
            if num_nest % 2 == 0:
                for next_nest in nests:
                    if next_nest == combined_nests:
                        continue
                    for comp_pair in combined_nests.competing_pairs:
                        if comp_pair in next_nest.nest:
                            competing_struc = True
                    if competing_struc == True:
                        combined_nests.length = len(combined_nests.nest)
                        next_nest.length = len(next_nest.nest)
                        if combined_nests.length >= next_nest.length:
                            for pair in next_nest.nest:
                                pk_pairs.append(pair)
                        else:
                            for pair in combined_nests.nest:
                                pk_pairs.append(pair)

    #print (pk_pairs, 'NESTED METHOD')
    return pk_pairs

def flag_1bp_prediction(pk_pairs):
    single_bp_prediction = []
    twoplus_prediction = []
    known_pairs = []
    if len(pk_pairs) == 1:
        single_bp_prediction = pk_pairs
        return [single_bp_prediction, []]
    for num, pair in enumerate(pk_pairs[:-1]):
        if (pair[0] != (pk_pairs[num+1][0]-1)) and (pair[0] != (pk_pairs[num+1][0]+1)) and (pair[1] != (pk_pairs[num+1][1]-1))  and (pair[1] != (pk_pairs[num+1][1]+1)):
            if pair not in known_pairs:
                single_bp_prediction.append(pair)
                known_pairs.append(pair)
        else:
            known_pairs.append(pair)
            known_pairs.append(pk_pairs[num+1])
    for pair in pk_pairs:
        if pair not in single_bp_prediction and pair not in twoplus_prediction:
            twoplus_prediction.append(pair)
    return [single_bp_prediction, twoplus_prediction]

def scorer_for_pks_check_1bp_predicitons(name, path_to_predicted_ct, path_to_ct='pseudobase.drct/known_pk_dot_and_ct/', spot_correction=1, spot_correction_predict=1):
    print(name,'five prime end of the region of interest')
    #os.chdir('pseudobase.drct')
    #print(os.getcwd())

    path_to_accepted_ct = path_to_ct+name+'.ct'
    #path_to_accepted_ct = path_to_ct
    accepted_file = open(path_to_accepted_ct, 'r')
    accepted_file_lines = accepted_file.readlines()
    accepted_pairs = []
    #the spot_correction variable serves to adjust for the format of the SPOT-RNA generated ct files, as they have an extra /n than the normal ct
    #this loop collects the pairs of i and j and adds them to a list
    for line in accepted_file_lines[spot_correction:]:
        i_and_j = []
        line_elements = line.split(' ')
        line_elements = ' '.join(line_elements).split()
        i_and_j.append(line_elements[0])
        i_and_j.append(line_elements[4])
        accepted_pairs.append(i_and_j)

    #the following block does the same as above but for the predicted file
    
    predicted_file = open(path_to_predicted_ct, 'r')
    predicted_file_lines = predicted_file.readlines()
    predicted_pairs = []
    for line in predicted_file_lines[spot_correction_predict:]:
        i_and_j = []
        line_elements = line.split(' ')
        line_elements = ' '.join(line_elements).split()
        i_and_j.append(line_elements[0])
        i_and_j.append(line_elements[4])
        predicted_pairs.append(i_and_j)

    sens_bp = 0
    sens_hits = 0
    ppv_bp = 0
    ppv_hits = 0
    sens_pk_bp = 0
    sens_pk_hits = 0
    ppv_pk_bp = 0
    ppv_pk_hits = 0


    sens_pk_bp_2 = 0
    sens_pk_hits_2 = 0
    ppv_pk_bp_2 = 0
    ppv_pk_hits_2 = 0

    #determines sensitivity for all base pairs, pseudoknots included
    num_true_pairs = 0
    num_agreement = 0
    already_paired = []
    already_counted = []
    for number, pair in enumerate(accepted_pairs):
        if int(pair[1]) != 0 and (pair[1]) not in already_counted:
            num_true_pairs += 1
            already_counted.append(pair[0])
        for compair in predicted_pairs:
            if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                if int(pair[1]) != 0:
                    num_agreement += 1
                    already_paired.append(pair[0])
                    already_paired.append(compair)

    sens_bp = num_true_pairs
    sens_hits = num_agreement
    try:
        print ('Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    except:
        pass
    #determines PPV for all base pairs, pseudoknots included
    num_true_pairs = 0
    num_agreement = 0
    already_paired = []
    already_counted = []
    for number, pair in enumerate(predicted_pairs):
        if int(pair[1]) != 0 and (pair[1]) not in already_counted:
            num_true_pairs += 1
            already_counted.append(pair[0])
        for compair in accepted_pairs:
            if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                if int(pair[1]) != 0:
                    num_agreement += 1
                    already_paired.append(pair[0])
                    already_paired.append(compair)
    ppv_bp = num_true_pairs
    ppv_hits = num_agreement
    if num_true_pairs != 0:
        print ('PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    else:
        print('No base pairs predicted')

    #this calls the pseudoknot base pair finder and builds lists of pairs
    #the following blocks repear a similar process as above to determine PPV and sensitivity. This should have been done by a function, but was not
    accepted_pk_pairs = nested_checker_1bp_prediciton(accepted_pairs)
    predicted_pk_pairs = nested_checker_1bp_prediciton(predicted_pairs)
    predicted_pk_pairs_flag = flag_1bp_prediction(predicted_pk_pairs)
    predicted_pk_pairs_1bp = predicted_pk_pairs_flag[0]
    predicted_pk_pairs_2plus_bp = predicted_pk_pairs_flag[1]

    if len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs_1bp) == 0:
        print ('No found pseudoknots in known and the predicted')
    elif len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs_1bp) != 0:
        num_agreement = 0
        num_true_pairs = len(predicted_pk_pairs_1bp)
        ppv_pk_bp = num_true_pairs
        print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    elif len(accepted_pk_pairs) != 0 and len(predicted_pk_pairs_1bp) == 0:
        num_agreement = 0
        num_true_pairs = len(accepted_pk_pairs)
        sens_pk_bp = num_true_pairs
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    else:
        #determines sensitivity for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []
        for number, pair in enumerate(accepted_pk_pairs):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                num_true_pairs += 1
                already_counted.append(pair[0])
            for num_com, compair in enumerate(predicted_pk_pairs_1bp):
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        sens_pk_bp = num_true_pairs
        sens_pk_hits = num_agreement
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))

        #determines the PPV for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []    
        for number, pair in enumerate(predicted_pk_pairs_1bp):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                #print(pair)
                num_true_pairs += 1
                already_counted.append(pair[0])
            for compair in accepted_pk_pairs:
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(pair)
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        ppv_pk_bp = num_true_pairs
        ppv_pk_hits = num_agreement
        if num_true_pairs != 0:
            print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
        else:
            print('No pseudoknots predicted')


    ### THIS IS FOR two+ bp redictions
    if len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs_2plus_bp) == 0:
        print ('No found pseudoknots in known and the predicted')
    elif len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs_2plus_bp) != 0:
        num_agreement = 0
        num_true_pairs = len(predicted_pk_pairs_2plus_bp)
        ppv_pk_bp_2 = num_true_pairs
        print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    elif len(accepted_pk_pairs) != 0 and len(predicted_pk_pairs_2plus_bp) == 0:
        num_agreement = 0
        num_true_pairs = len(accepted_pk_pairs)
        sens_pk_bp_2 = num_true_pairs
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    else:
        #determines sensitivity for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []
        for number, pair in enumerate(accepted_pk_pairs):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                num_true_pairs += 1
                already_counted.append(pair[0])
            for num_com, compair in enumerate(predicted_pk_pairs_2plus_bp):
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        sens_pk_bp_2 = num_true_pairs
        sens_pk_hits_2 = num_agreement
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))

        #determines the PPV for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []    
        for number, pair in enumerate(predicted_pk_pairs_2plus_bp):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                #print(pair)
                num_true_pairs += 1
                already_counted.append(pair[0])
            for compair in accepted_pk_pairs:
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(pair)
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        ppv_pk_bp_2 = num_true_pairs
        ppv_pk_hits_2 = num_agreement
        if num_true_pairs != 0:
            print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
        else:
            print('No pseudoknots predicted')


    return [[sens_hits, sens_bp, ppv_hits, ppv_bp, sens_pk_hits, sens_pk_bp, ppv_pk_hits, ppv_pk_bp], [sens_hits, sens_bp, ppv_hits, ppv_bp, sens_pk_hits_2, sens_pk_bp_2, ppv_pk_hits_2, ppv_pk_bp_2]]


#This function serves to score pseudoknot base pairs in a given ct file.
#it compares by default the predicted ct file with the known ct file from pseudobase, though the known structure can be changed.
def scorer_for_pks(name, path_to_predicted_ct, path_to_ct='pseudobase.drct/known_pk_dot_and_ct/', spot_correction=1, spot_correction_predict=1):
    print(name,'five prime end of the region of interest')
    #os.chdir('pseudobase.drct')
    #print(os.getcwd())

    #path_to_accepted_ct = path_to_ct+name+'.ct'
    path_to_accepted_ct = path_to_ct
    accepted_file = open(path_to_accepted_ct, 'r')
    accepted_file_lines = accepted_file.readlines()
    accepted_pairs = []
    #the spot_correction variable serves to adjust for the format of the SPOT-RNA generated ct files, as they have an extra /n than the normal ct
    #this loop collects the pairs of i and j and adds them to a list
    for line in accepted_file_lines[spot_correction:]:
        i_and_j = []
        line_elements = line.split(' ')
        line_elements = ' '.join(line_elements).split()
        i_and_j.append(line_elements[0])
        i_and_j.append(line_elements[4])
        accepted_pairs.append(i_and_j)

    #the following block does the same as above but for the predicted file
    predicted_file = open(path_to_predicted_ct, 'r')
    predicted_file_lines = predicted_file.readlines()
    predicted_pairs = []
    for line in predicted_file_lines[spot_correction_predict:]:
        i_and_j = []
        line_elements = line.split(' ')
        line_elements = ' '.join(line_elements).split()
        i_and_j.append(line_elements[0])
        i_and_j.append(line_elements[4])
        predicted_pairs.append(i_and_j)

    sens_bp = 0
    sens_hits = 0
    ppv_bp = 0
    ppv_hits = 0
    sens_pk_bp = 0
    sens_pk_hits = 0
    ppv_pk_bp = 0
    ppv_pk_hits = 0

    #determines sensitivity for all base pairs, pseudoknots included
    #number of actual true base pairs there are in the accepted structure
    num_true_pairs = 0
    #number of true base pairs found in the predicted structure
    num_agreement = 0
    already_paired = []
    already_counted = []
    for number, pair in enumerate(accepted_pairs):
        #first ensures that it is a base pair, the second term ensures the i wasn't already counted for a j
        if int(pair[1]) != 0 and (pair[1]) not in already_counted:
            num_true_pairs += 1
            already_counted.append(pair[0])
        for compair in predicted_pairs:
            if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                if int(pair[1]) != 0:
                    num_agreement += 1
                    already_paired.append(pair[0])
                    already_paired.append(compair)

    sens_bp = num_true_pairs
    sens_hits = num_agreement
    try:
        print ('Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    except:
        pass
    #determines PPV for all base pairs, pseudoknots included
    num_true_pairs = 0
    num_agreement = 0
    already_paired = []
    already_counted = []
    for number, pair in enumerate(predicted_pairs):
        if int(pair[1]) != 0 and (pair[1]) not in already_counted:
            num_true_pairs += 1
            already_counted.append(pair[0])
        for compair in accepted_pairs:
            if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                if int(pair[1]) != 0:
                    num_agreement += 1
                    already_paired.append(pair[0])
                    already_paired.append(compair)
    ppv_bp = num_true_pairs
    ppv_hits = num_agreement
    if num_true_pairs != 0:
        print ('PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    else:
        print('No base pairs predicted')

    #this calls the pseudoknot base pair finder and builds lists of pairs
    #the following blocks repear a similar process as above to determine PPV and sensitivity. This should have been done by a function, but was not
    accepted_pk_pairs = nested_checker(accepted_pairs)
    predicted_pk_pairs = nested_checker(predicted_pairs)

    if len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs) == 0:
        print ('No found pseudoknots in known and the predicted')
    elif len(accepted_pk_pairs) == 0 and len(predicted_pk_pairs) != 0:
        num_agreement = 0
        num_true_pairs = len(predicted_pk_pairs)
        ppv_pk_bp = num_true_pairs
        print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    elif len(accepted_pk_pairs) != 0 and len(predicted_pk_pairs) == 0:
        num_agreement = 0
        num_true_pairs = len(accepted_pk_pairs)
        sens_pk_bp = num_true_pairs
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
    else:
        #determines sensitivity for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []
        for number, pair in enumerate(accepted_pk_pairs):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                num_true_pairs += 1
                already_counted.append(pair[0])
            for num_com, compair in enumerate(predicted_pk_pairs):
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        sens_pk_bp = num_true_pairs
        sens_pk_hits = num_agreement
        print ('Pseudoknot Base Pair Sensitivity: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))

        #determines the PPV for pseudoknot base pairs
        num_true_pairs = 0
        num_agreement = 0
        already_paired = []
        already_counted = []    
        for number, pair in enumerate(predicted_pk_pairs):
            if int(pair[1]) != 0 and (pair[1]) not in already_counted:
                #print(pair)
                num_true_pairs += 1
                already_counted.append(pair[0])
            for compair in accepted_pk_pairs:
                if ((pair[1] not in already_paired) and (compair not in already_paired) and ( (int(pair[0]) == int(compair[0]) and (int(compair[1]) in range(int(pair[1])-1, int(pair[1])+2))) or (int(pair[1] == int(compair[1])) and (int(compair[0]) in range(int(pair[0])-1, int(pair[0])+2)) ) )):    
                    if int(pair[1]) != 0:
                        #print(pair)
                        #print(compair, '\n')
                        num_agreement += 1
                        already_paired.append(pair[0])
                        already_paired.append(compair)
        #num_agreement /= 2
        #num_true_pairs /= 2
        ppv_pk_bp = num_true_pairs
        ppv_pk_hits = num_agreement
        if num_true_pairs != 0:
            print ('Pseudoknot Base Pair PPV: %d / %d = %s' % (num_agreement, num_true_pairs, str(float(num_agreement)/float(num_true_pairs)*100)[:5]+'%'))
        else:
            print('No pseudoknots predicted')


    return [sens_hits, sens_bp, ppv_hits, ppv_bp, sens_pk_hits, sens_pk_bp, ppv_pk_hits, ppv_pk_bp]

#this quick function determines if two base pairs are nested or not
def is_nested(nest, new_pair):
    new_duo = [0,0]
    if new_pair[1] < new_pair[0]:
        new_duo[0] = new_pair[1]
        new_duo[1] = new_pair[0]
    else:
        new_duo = new_pair
    if new_duo[0] > nest[0] and new_duo[1] < nest[1]:
        return True
    if (new_duo[0] > nest[0] and new_duo[1] > nest[1]) or (new_duo[0] < nest[0] and new_duo[1] < nest[1]):
        return False

#this function determines if two pairs are nested and continuous without bulges or loops
def is_continuous_nest(nest):
    continuousR = 0
    continuousL = 0
    for number, pair in enumerate(nest[2]):
        if number == len(nest[2]) -1:
            break
        if pair[0] + 1 == nest[2][number+1]:
            continuousR = 1
        if pair[1] - 1 == nest[2][number+1]:
            continuousL = 1
    return continuousR + continuousL

#unclear why another function that does the exact thing as the previous one was needed
def is_continuous_pknest(pk_nest):
    continuousR = 0
    continuousL = 0
    for number, pair in enumerate(pk_nest[0]):
        if number < len(pk_nest[0])-1:
            if pair[0] + 1 == pk_nest[0][number+1]:
                continuousR = 1
            if pair[1] - 1 == pk_nest[0][number+1]:
                continuousL = 1
    return continuousR + continuousL    

#the function tests the nested_checker function using the probknot output cts
def nested_checker_checker():
    os.chdir('pseudobase.drct/probknot_massive')
    for file in os.listdir():
        if '.ct' in file:
            print(file)
            accepted_file = open(file, 'r')
            accepted_file_lines = accepted_file.readlines()
            accepted_pairs = []
            for line in accepted_file_lines[1:]:
                i_and_j = []
                line_elements = line.split(' ')
                line_elements = ' '.join(line_elements).split()
                i_and_j.append(line_elements[0])
                i_and_j.append(line_elements[4])
                accepted_pairs.append(i_and_j)
            nested_checker(accepted_pairs)

class nest_class:
    def __init__(self):
        self.nest = []
        self.length = 0
        self.competing_pairs = []

class pair_class:
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.closed = False
        self.pair = [i, j]

#this function is part of the pseudoknot-specific scoring method
#this finds the non-nested base pairs that are known as pseudoknots
def nested_checker(accepted_pairs):
    nests = []
    pk_pairs = []
    known_bases = []
    end = False
    #this creates continuous 'nests' out of the accepted pairs
    for line, pair in enumerate(accepted_pairs):
        pair[0] = int(pair[0])
        pair[1] = int(pair[1])
        if pair[1] != 0 and pair[0] not in known_bases:
            #starts the first nest
            if nests == []:
                temp_nest = nest_class()
                temp_nest.nest.append(pair)
                nests.append(temp_nest)
            #1st term: is base adjacent to the preveious one, 2nd: is the j base after the previous j
            #this determines if they're "nested base pairs"
            elif pair[0] - 1 == accepted_pairs[line-1][0] and pair[1]+1 == accepted_pairs[line-1][1]:
                temp_nest.nest.append(pair)
            #this starts a new class for non-nests
            else:
                temp_nest = nest_class()
                temp_nest.nest.append(pair)
                nests.append(temp_nest)
            known_bases.append(pair[1])

    nests_to_remove = []

    #this removes nests that are totally self-nested (definition of a hairpin), and therefore do not participate in pseudoknots
    for nest in nests:
        remove = True
        for pair in accepted_pairs:
            #selects a pair with i after the  nest's final i term
            if pair[0] > nest.nest[-1][0]:
                #this conditional is the key... if the pair is outside of the nest and does not form a pk, then the nest is non-pseudoknotted
                if pair[1] != 0 and pair[0] <= nest.nest[-1][1]:
                    if is_nested(nest.nest[-1], pair) == False:
                        added_pair = [0,0]
                        if pair[1] < pair[0]:
                            added_pair[0] = pair[1]
                            added_pair[1] = pair[0]
                        else:
                            added_pair = pair
                        nest.competing_pairs.append(added_pair)
                        remove = False
        if remove == True:
            nests_to_remove.append(nest)
    for nest in nests_to_remove:
        nests.remove(nest)

    #this serves to combine nests that are separated by bulges or breaks and are nested
    for num_ses, suspect_nest in enumerate(nests):
        for nest_2 in nests[num_ses:]:
            same_comp = False
            if is_nested(suspect_nest.nest[-1], nest_2.nest[-1]) == True:
                for comp_base in suspect_nest.competing_pairs:
                    #ensures detected pseudoknots are kept
                    if comp_base in nest_2.competing_pairs:
                        same_comp = True
                if same_comp == True:
                    for pair_base in nest_2.nest:
                        suspect_nest.nest.append(pair_base)
                    nests.remove(nest_2)
    
    #this serves to combine two nests that are borken by the same pseudoknot
    if len(nests) % 2 != 0:
        for num_nest, nest in enumerate(nests):
            add_pairs = False
            for next_nest in nests[num_nest+1:]:
                for next_pair in next_nest.competing_pairs:
                    if next_pair in nest.competing_pairs: 
                        add_pairs = True
                        break
                if add_pairs == True:
                    for nest_next_pair in next_nest.nest:
                        nest.nest.append(nest_next_pair)
                    nests.remove(next_nest)

    #this is the final determiner of what combined nests of stuctures are pseudoknots and what are nested
    if len(nests) != 0 and len(nests) != 1:
        for num_nest, combined_nests in enumerate(nests):
            competing_struc = False
            if num_nest % 2 == 0:
                for next_nest in nests:
                    if next_nest == combined_nests:
                        continue
                    for comp_pair in combined_nests.competing_pairs:
                        if comp_pair in next_nest.nest:
                            competing_struc = True
                    if competing_struc == True:
                        combined_nests.length = len(combined_nests.nest)
                        next_nest.length = len(next_nest.nest)
                        if combined_nests.length >= next_nest.length:
                            for pair in next_nest.nest:
                                pk_pairs.append(pair)
                        else:
                            for pair in combined_nests.nest:
                                pk_pairs.append(pair)

    #print (pk_pairs, 'NESTED METHOD')
    return pk_pairs

#OBSOLETE... this function previously foudn pseudoknot basepairs. Replaced by nested_checker
def pseudoknot_base_pairs(accepted_pairs):
    nest_1 = [[0,0], 0, []]
    nest_2 = [[], 0]
    nest_3 = [[0,0], 0]
    pk_nest = [[], 0]
    running_nest = [0,0]
    known_bases = []
    pk_pairs = []
    end = False
    #this loop runs through all the known pairs
    for pair in accepted_pairs:
        is_pk_the_pk = True
        #if pk_nest[1] <= nest_1[1] and is_continuous_pknest(pk_nest) <= is_continuous_nest(nest_1):
        #    temp_nest = nest_1
         #   nest_1 = [[0,0], 0, []]
        #    for pk_thing in pk_nest[2:]:
        #        nest_1[2].append(pk_thing)
        #        nest_1[1] = pk_nest[1]
        #        nest_1[0] = pk_nest[0]
        #    pk_nest[1] = temp_nest[1]
        #    for i in temp_nest[2]:
        #        pk_nest[0].append(i)
        #print(pair)
        #print(nest_1, 'AAAAAA')
        #print(running_nest, 'BBBBB')
        #print(pk_nest, 'CCCCC')
        pair[0] = int(pair[0])
        pair[1] = int(pair[1])
        #print(nest_1)
        #this conditionals checks to see if the loop has reached the end of the sequence
        if pair == accepted_pairs[-1]:
            end = True
        #this conditional ensures that the pair being analyzed has not been already analyzed
        if pair[0] not in known_bases or end == True:
            #if nest_1[0] == [0,0] and pair[1] != 0:
            #if it is the first base pair
            if nest_1[0] == [0,0] and pair[1] != 0:
                nest_1[0] = pair
                nest_1[1] += 1
                nest_1[2].append(pair)
                running_nest = pair
            #elif (pair[0] > running_nest[0]) and (pair[1] < running_nest[1]) and pair[1] != 0 and end == False:
            #if the new pair is nested in the old pair
            elif is_nested(running_nest, pair) == True and pair[1] != 0 and end == False:
                running_nest = pair
                nest_1[1] += 1
                nest_1[2].append(pair)
            #if the new pair is outside the nest of the old pair completely
            elif pair[0] > nest_1[0][1] or end == True:
                #if the pk pairs are less in length to the real nest or less contiguous
                if pk_nest[1] <= nest_1[1] and is_continuous_pknest(pk_nest) <= is_continuous_nest(nest_1):
                    for pk_bp in pk_nest[0]:
                        pk_pairs.append(pk_bp)
                    pk_nest = [[], 0]
                elif pk_nest[1]-1 == nest_1[1] and is_continuous_pknest(pk_nest) <= is_continuous_nest(nest_1):
                    for pk_bp in pk_nest[0]:
                        pk_pairs.append(pk_bp)
                    pk_nest = [[], 0]
                else:

                    #pk_pairs.append(nest_1[2])
                    for pk_bp in nest_1[2]:
                        if pk_bp not in pk_pairs:
                            pk_pairs.append(pk_bp)
                    is_pk_the_pk = False
                #conditional to end the accepted pairs list
                if pair[1] == 0 and end == True:
                    nest_1[0] = [[0,0], 0, []] 
                    running_nest = [0,0]
                    pk_nest = [[], 0]    
                else:
                    if is_pk_the_pk == True:
                        nest_1[0] = pair
                        running_nest = pair
                        nest_1[1] = 1
                        nest_1[2] = [pair]
                        pk_nest = [[], 0]
                    else:
                        running_nest = nest_1[2][-1]
                        #nest_1 = [[0,0], 0, []]
                        if is_nested(running_nest, pair) == False:
                            nest_1[0] = pk_nest[0][0]
                            nest_1[1] = pk_nest[1]
                            nest_1[2] = []
                            for bases in pk_nest[0]:
                                nest_1[2].append(bases)
                            pk_nest[0] = []
                            pk_nest[0].append(pair)
                            pk_nest[1] = 1
            elif is_nested(running_nest, pair) == False and pair[1] != 0 and end == False:
                if is_nested(nest_1[0], pair) == False:
                    #print(nest_1[0])
                    pk_nest[0].append(pair)
                    pk_nest[1] += 1
                elif pair[0] > running_nest[1]:
                    if nest_1[1] >= pk_nest[1]:
                        nest_1[1] += 1
                        nest_1[2].append(pair)
                        #print(nest_1)
                    else:
                        pk_nest[1] += 1
                        pk_nest[0].append(pair)
                else:
                    pk_nest[0].append(pair)
                    pk_nest[1] += 1
        known_bases.append(pair[0])
        known_bases.append(pair[1])

    #flat_pk_pairs = []
    #for sublist in pk_pairs:
    #     for item in sublist:
    #         flat_pk_pairs.append(item)
    #print (pk_pairs)
    print (pk_pairs, 'OLD METHOD')
    return pk_pairs

#this function makes fasta files
def make_fasta_files(name, sequence):
    temp_file_name = name+'.fasta'
    temp_file = open(temp_file_name, 'w')
    temp_file.write('>'+name+'\n'+sequence)
    temp_file.close()
    return temp_file_name

#this function ran ipknot on all the postive control sequences and compared the results to the known database
def massive_test_ipknot():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('ipknot_massive')
    except:
        pass
    os.chdir('ipknot_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        fasta_file_name = make_fasta_files(name, sequence)
        ipknot_output_txt_name = name +'.txt'
        ipknot_output_ct_name = name +'.ct'
        subprocess.run('ipknot %s > %s' % (fasta_file_name, ipknot_output_txt_name), shell=True)
        ipknot_output_txt = open(ipknot_output_txt_name, 'r')
        output_data = ipknot_output_txt.readlines()
        new_string = ''
        for line in output_data[1:]:
            new_string += line
        ipknot_output_txt.close()
        ipknot_output_txt = open(ipknot_output_txt_name, 'w')
        ipknot_output_txt.write(new_string)
        ipknot_output_txt.close()
        subprocess.run('dot2ct %s %s' % (ipknot_output_txt_name, ipknot_output_ct_name), shell=True)
        comparison = scorer_for_pks(name, ipknot_output_ct_name)
        #number = export_to_excel('IPKnot', name, comparison, number)
        for place, value in enumerate(comparison):
            final_values[place] += value
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values[0], final_values[1], str(float(final_values[0])/float(final_values[1])*100)+'%'))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values[2], final_values[3], str(float(final_values[2])/float(final_values[3])*100)+'%'))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values[4], final_values[5], str(float(final_values[4])/float(final_values[5])*100)+'%'))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values[6], final_values[7], str(float(final_values[6])/float(final_values[7])*100)+'%'))
    except:
        print('No pseudoknot base pairs predicted')

#this function ran probknot on all the postive control sequences and compared the results to the known database
def massive_test_prob_knot():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('probknot_massive')
    except:
        pass
    os.chdir('probknot_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        fasta_file_name = make_fasta_files(name, sequence)
        probknot_output_ct = name +'.ct'
        subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, probknot_output_ct), shell=True)
        comparison = scorer_for_pks(name, probknot_output_ct)
        #number = export_to_excel('ProbKnot', name, comparison, number)
        for place, value in enumerate(comparison):
            final_values[place] += value
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values[0], final_values[1], str(float(final_values[0])/float(final_values[1]))))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values[2], final_values[3], str(float(final_values[2])/float(final_values[3]))))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values[4], final_values[5], str(float(final_values[4])/float(final_values[5]))))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values[6], final_values[7], str(float(final_values[6])/float(final_values[7]))))
    except:
        pass

#this function ran spot-rna on all the postive control sequences and compared the results to the known database
def massive_test_spot_rna():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('spotrna_massive')
    except:
        pass
    os.chdir('spotrna_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        try:
            os.mkdir(name+'.drct')
        except:
            pass
        os.chdir(name+'.drct')
        fasta_file_name = make_fasta_files(name, sequence)
        spot_rna_output_ct = name +'.ct'
        while os.getcwd() != '/':
                os.chdir('../')
        os.chdir('home/pfors/programs/SPOT-RNA')
        spot_rna_output_path= '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Reseach Lab/ScanFold/ScanFold_working_directory/pseudobase.drct/spotrna_massive/'+name+'.drct/'
        spot_rna_output_path_huh= '/mnt/c/Users/pfors/Desktop/Research\ Lab/Bevilacqua\ Reseach\ Lab/ScanFold/ScanFold_working_directory/pseudobase.drct/spotrna_massive/'+name+'.drct/'
        #subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+name+'.fasta --outputs '+spot_rna_output_path_huh+' --plots True', shell=True)
        os.chdir(spot_rna_output_path)
        os.chdir('../')
        #print(os.getcwd())
        comparison = scorer_for_pks(name, name+'.drct/'+spot_rna_output_ct ,'../known_pk_dot_and_ct/',1,2)
        #number = export_to_excel('SPOT-RNA', name, comparison, number)
        #subprocess.run('ls')
        #s.chdir('../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/")
        #'../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/"+str(pk.five_bond)+'.ct'
        #is_spot_sig = SPOT_RNA_significance('../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/"+str(pk.five_bond)+'.ct', str(pk.five_prime_boundary)+'.ct', spot_rna_output_name, degree_support)
        #pk.num_pks = num_pks
        #return is_spot_sig
        for place, value in enumerate(comparison):
            final_values[place] += value
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values[0], final_values[1], str(float(final_values[0])/float(final_values[1])*100)+'%'))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values[2], final_values[3], str(float(final_values[2])/float(final_values[3])*100)+'%'))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values[4], final_values[5], str(float(final_values[4])/float(final_values[5])*100)+'%'))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values[6], final_values[7], str(float(final_values[6])/float(final_values[7])*100)+'%'))
    except:
        pass

def massive_test_spot_rna_1bp():
    os.chdir('pseudobase.drct')
    try:
        os.mkdir('spotrna_massive')
    except:
        pass
    os.chdir('spotrna_massive')
    names_file_name = '../stored_pseudobase_extractions_ALL.txt'
    names_file = open(names_file_name, 'r')
    final_values_1bp = [0,0,0,0,0,0,0,0]
    final_values_2bp = [0,0,0,0,0,0,0,0]
    names_lines = names_file.readlines()
    number = 0
    for line in names_lines[1:]:
        elements = line.split('\t')
        name = elements[0]
        sequence = elements[8]
        if '.' in sequence or 'Y' in sequence:
            continue
        try:
            os.mkdir(name+'.drct')
        except:
            pass
        os.chdir(name+'.drct')
        fasta_file_name = make_fasta_files(name, sequence)
        spot_rna_output_ct = name +'.ct'
        while os.getcwd() != '/':
                os.chdir('../')
        os.chdir('home/pfors/programs/SPOT-RNA')
        spot_rna_output_path= '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Reseach Lab/ScanFold/ScanFold_working_directory/pseudobase.drct/spotrna_massive/'+name+'.drct/'
        spot_rna_output_path_huh= '/mnt/c/Users/pfors/Desktop/Research\ Lab/Bevilacqua\ Reseach\ Lab/ScanFold/ScanFold_working_directory/pseudobase.drct/spotrna_massive/'+name+'.drct/'
        #subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+name+'.fasta --outputs '+spot_rna_output_path_huh+' --plots True', shell=True)
        os.chdir(spot_rna_output_path)
        os.chdir('../')
        #print(os.getcwd())
        comparison = scorer_for_pks_check_1bp_predicitons(name, name+'.drct/'+spot_rna_output_ct ,'../known_pk_dot_and_ct/',1,2)
        #number = export_to_excel('SPOT-RNA', name, comparison, number)
        #subprocess.run('ls')
        #s.chdir('../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/")
        #'../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/"+str(pk.five_bond)+'.ct'
        #is_spot_sig = SPOT_RNA_significance('../../Umb_d/NC_045512.2_SARS_CoV_2_ref_seq.drct/NC_0455122Severeacut.drct/SPOT_RNA_files.drct/'+str(pk.five_bond)+".drct/"+str(pk.five_bond)+'.ct', str(pk.five_prime_boundary)+'.ct', spot_rna_output_name, degree_support)
        #pk.num_pks = num_pks
        #return is_spot_sig
        for place, value in enumerate(comparison[0]):
            final_values_1bp[place] += value
        for place, value in enumerate(comparison[1]):
            final_values_2bp[place] += value    
    print('\n\n\n\n')
    print('FOR 1bp PREDICTIONS')
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values_1bp[0], final_values_1bp[1], str(float(final_values_1bp[0])/float(final_values_1bp[1])*100)+'%'))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values_1bp[2], final_values_1bp[3], str(float(final_values_1bp[2])/float(final_values_1bp[3])*100)+'%'))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values_1bp[4], final_values_1bp[5], str(float(final_values_1bp[4])/float(final_values_1bp[5])*100)+'%'))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values_1bp[6], final_values_1bp[7], str(float(final_values_1bp[6])/float(final_values_1bp[7])*100)+'%'))
    except:
        pass
    print('\n\n\n\n')
    print('FOR 2+bp PREDICTIONS')
    print ('All Base Pair Sensitivity = %d / %d = %s' % (final_values_2bp[0], final_values_2bp[1], str(float(final_values_2bp[0])/float(final_values_2bp[1])*100)+'%'))
    print ('All Base Pair PPV = %d / %d = %s' % (final_values_2bp[2], final_values_2bp[3], str(float(final_values_2bp[2])/float(final_values_2bp[3])*100)+'%'))
    print ('All Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (final_values_2bp[4], final_values_2bp[5], str(float(final_values_2bp[4])/float(final_values_2bp[5])*100)+'%'))
    try:
        print ('All Pseudoknot Base Pair PPV = %d / %d = %s' % (final_values_2bp[6], final_values_2bp[7], str(float(final_values_2bp[6])/float(final_values_2bp[7])*100)+'%'))
    except:
        pass
    print('\n\n\n\n')




'''
#this program reported the overlapping hits of structured pseudoknots.
#OBSOLETE
#determined to be defunct and error-prone
def export_to_excel(test_method, name, comparison_list, number):
    excel_name = 'Overlap_of_pseudoknot_programs_spotrna_probknot_ipknot.xlsx'
    #excel_book = Workbook(excel_name)
    #excel_book.save(excel_name)
    excel_book=load_workbook(excel_name)
    sheet = excel_book.active
    test = 1
    ppv = 0
    if test_method == 'SPOT-RNA':
        test = 2
    elif test_method == 'ProbKnot':
        test = 3
    elif test_method == 'IPKnot':
        test = 4
    if comparison_list[7]  != 0:
        ppv = 1

    sheet.cell(row=1, column=1).value = 'Name'
    sheet.cell(row=1, column=2).value = 'SPOT-RNA'
    sheet.cell(row=1, column=3).value = 'ProbKnot'
    sheet.cell(row=1, column=4).value = 'IPKnot'

    sheet.cell(row=number+2, column=1).value = name
    sheet.cell(row=number+2, column=test).value = ppv
    excel_book.save(excel_name)

    number += 1
    return number

#this program reported the overlapping hits of structured pseudoknots.
#OBSOLETE
#determined to be defunct and error-prone
def read_excel_overlap():
    excel_book_name =  'Overlap_of_pseudoknot_programs_spotrna_probknot_ipknot.xlsx'
    excel_book = load_workbook(excel_book_name)
    sheet = excel_book.active

    total_spot = 0
    spot_and_prob = 0
    spot_and_ip = 0
    spot_and_all = 0

    total_prob = 0
    prob_and_ip = 0

    total_ip = 0
    
    row_num = 2
    while row_num < 370:
        if sheet.cell(row=row_num, column=1).value == None:
            print('\n +end list +\n')
            break
        print(row_num)
        total_spot += sheet.cell(row=row_num,column=2).value
        total_prob += sheet.cell(row=row_num,column=3).value
        total_ip += sheet.cell(row=row_num,column=4).value
        if sheet.cell(row=row_num,column=2).value == 1 and sheet.cell(row=row_num,column=3).value == 1:
            spot_and_prob += 1
        if sheet.cell(row=row_num,column=2).value == 1 and sheet.cell(row=row_num,column=4).value == 1:
            spot_and_ip += 1
        if sheet.cell(row=row_num,column=2).value == 1 and sheet.cell(row=row_num,column=3).value == 1 and sheet.cell(row=row_num,column=4).value == 1:
            spot_and_all += 1
        if sheet.cell(row=row_num,column=3).value == 1 and sheet.cell(row=row_num,column=4).value == 1:        
            prob_and_ip += 1
        row_num += 1

    print('total spot', total_spot, ' | total prob',total_prob,' | total ip',total_ip)
    print('SPOT and Prob overlap', spot_and_prob, ' | SPOT and ip overlap', spot_and_ip, ' | Spot and all overlap', spot_and_all, ' | prob and ip overlap', prob_and_ip)
    spot_in_prob = float(spot_and_prob) / float(total_spot)
    print('spot in prob', spot_in_prob)
    prob_in_spot = float(spot_and_prob) / float(total_prob)
    print('prob_in_spot',prob_in_spot)

    spot_in_ip = float(spot_and_ip) / float(total_spot)
    print('spot_in_ip',spot_in_ip)
    ip_in_spot = float(spot_and_ip) / float(total_ip)
    print('ip_in_spot',ip_in_spot)

    prob_in_ip = float(prob_and_ip) / float(total_prob)
    print('prob_in_ip',prob_in_ip)
    ip_in_prob = float(prob_and_ip) / float(total_ip)
    print('ip_in_prob',ip_in_prob)

    spot_in_all = float(spot_and_all) / float(total_spot)
    print('spot_in_all',spot_in_all)
    prob_in_all = float(spot_and_all) / float(total_prob)
    print('prob_in_all',prob_in_all)
    ip_in_all = float(spot_and_all) / float(total_ip)
    print('ip_in_all',ip_in_all)
'''
#this function calls all three pseudoknot programs and collects their outputs
def all_pk_programs(fasta_file_name):
    fasta_file = open(fasta_file_name, 'r')
    lines = fasta_file.readlines()
    seq_name = lines[0]
    seq_name = seq_name.replace(' ', '_')
    seq_name = seq_name.replace('\n', '')
    seq_name = seq_name.replace('>', '')
    print(seq_name)
    fasta_file.close()
    probknot_output_name = seq_name +'_probknot.ct'
    subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, probknot_output_name), shell=True)

    ipknot_output_txt_name = seq_name +'_ipknot.txt'
    ipknot_output_ct_name = seq_name +'_ipknot.ct'
    subprocess.run('ipknot %s > %s' % (fasta_file_name, ipknot_output_txt_name), shell=True)
    ipknot_output_txt = open(ipknot_output_txt_name, 'r')
    output_data = ipknot_output_txt.readlines()
    new_string = ''
    for line in output_data[1:]:
        new_string += line
    ipknot_output_txt.close()
    ipknot_output_txt = open(ipknot_output_txt_name, 'w')
    ipknot_output_txt.write(new_string)
    ipknot_output_txt.close()
    subprocess.run('dot2ct %s %s' % (ipknot_output_txt_name, ipknot_output_ct_name), shell=True)

    spot_rna_output_ct = seq_name +'_spotrna.ct'
    while os.getcwd() != '/':
            os.chdir('../')
    os.chdir('home/pfors/programs/SPOT-RNA')
    spot_rna_output_path= '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Reseach Lab/ScanFold/ScanFold_working_directory/'
    spot_rna_output_path_huh= '/mnt/c/Users/pfors/Desktop/Research\ Lab/Bevilacqua\ Reseach\ Lab/ScanFold/ScanFold_working_directory/'
    subprocess.run('python3 SPOT-RNA.py --inputs '+spot_rna_output_path_huh+fasta_file_name+' --outputs '+spot_rna_output_path_huh+' --plots True', shell=True)
    os.chdir(spot_rna_output_path)
    os.chdir('../')

#this function served to compare SNP pseudoknot prediction to the accepted predicted ct structure to determine if SPOT-RNA was capable of finding pseudoSNitches
def pseudoSNitch_analysis_SPOT(fasta_file_name):
    comparison_output_file_name = fasta_file_name +'_comparison_output.txt'
    comparison_output_file = open(comparison_output_file_name, 'a')
    bases_of_interest = [[3, 'G'], [11, 'A'],[12, 'C'],[15, 'U'],[16, 'G'],[17, 'C'],[19, 'A' ],[28 ,'C'],[29 ,'C'],[30 ,'C'],[33, 'G'],[34, 'C'],[40, 'U']]
    #base_format = [[base_number, identity, expected_outcome]]
    #expected_outcome = 'loss of pk', 'no change', 'change in secondary structure'
    fasta_file = open(fasta_file_name, 'r')
    fasta_lines = fasta_file.readlines()
    sequence = fasta_lines[1]
    name = fasta_lines[0][1:].replace('\n', '')
    while os.getcwd() != '/':
        os.chdir('../')
    os.chdir('home/pfors/programs/SPOT-RNA')
    subprocess.run('python3 SPOT-RNA.py --inputs %s --outputs outputs/' % (fasta_file_name),shell=True)
    print(name)
    for base in bases_of_interest:
        print(base)
        four_bases = ['A', 'U', 'C', 'G']
        for num_char, char in enumerate(sequence):
            if char == base[1] and num_char+1 == base[0]:
                print ('FOUND')
                print(char, num_char)
                four_bases.remove(char)
                for mut in four_bases:
                    mut_sequence = sequence[:base[0]-1]+mut+sequence[base[0]:]
                    mut_name = base[1]+str(base[0])+mut
                    total_mut_name = name+mut_name
                    print(total_mut_name)
                    temp_fasta = make_fasta_files(total_mut_name, mut_sequence)
                    temp_ct_output_name = total_mut_name+'_probknot_output.ct'
                    subprocess.run('python3 SPOT-RNA.py --inputs %s --outputs outputs/' % (temp_fasta),shell=True)
                    comparison = scorer_for_pks(name, 'outputs/'+total_mut_name+'.ct', 'outputs/',2,2)
                    comparison_string = ('All Base Pair Sensitivity = %d / %d = %s' % (comparison[0], comparison[1], str(float(comparison[0])/float(comparison[1])*100)+'%'))+('\nAll Base Pair PPV = %d / %d = %s' % (comparison[2], comparison[3], str(float(comparison[2])/float(comparison[3])*100)+'%'))+('\nAll Pseudoknot Base Pair Sensitivity = %d / %d = %s' % (comparison[4], comparison[5], str(float(comparison[4])/float(comparison[5])*100)+'%'))
                    try:
                        comparison_string += ('\nAll Pseudoknot Base Pair PPV = %d / %d = %s' % (comparison[6], comparison[7], str(float(comparison[6])/float(comparison[7])*100)+'%'))
                    except:
                        pass
                    comparison_output_file.write(total_mut_name+'\n'+comparison_string+'\n'+'\n')
            else:
                print ('Invalid base', num_char, char)

#this function served to compare SNP pseudoknot prediction to the accepted prediction of probknot to determine its ability to predict pseudoSNtiches
def pseudoSNitch_analysis_probknot(fasta_file_name):
    bases_of_interest = [[4, 'C']]
    #base_format = [[base_number, identity, expected_outcome]]
    #expected_outcome = 'loss of pk', 'no change', 'change in secondary structure'
    fasta_file = open(fasta_file_name, 'r')
    fasta_lines = fasta_file.readlines()
    sequence = fasta_lines[1]
    name = fasta_lines[0][1:].replace('\n', '')
    subprocess.run('ProbKnot --sequence %s %s' % (fasta_file_name, name+'.ct'),shell=True)
    print(name)
    for base in bases_of_interest:
        print(base)
        four_bases = ['A', 'U', 'C', 'G']
        for num_char, char in enumerate(sequence):
            if char == base[1] and num_char+1 == base[0]:
                print ('FOUND')
                print(char, num_char)
                four_bases.remove(char)
                for mut in four_bases:
                    mut_sequence = sequence[:base[0]-1]+mut+sequence[base[0]:]
                    mut_name = base[1]+str(base[0])+mut
                    total_mut_name = name+mut_name
                    print(total_mut_name)
                    temp_fasta = make_fasta_files(total_mut_name, mut_sequence)
                    temp_ct_output_name = total_mut_name+'_probknot_output.ct'
                    subprocess.run('ProbKnot --sequence %s %s' % (temp_fasta, temp_ct_output_name),shell=True)
                    comparison = scorer_for_pks(name, temp_ct_output_name, '')

                    #print(mut_sequence)
    #print(mutant_sequences)


#extract_pseudobase(html)
#extract_all('http://www.ekevanbatenburg.nl/PKBASE/PKBGETCLS.HTML')
#pseudobase_2_ct()
#clean_dir()
#check_for_doubles()
#scorer_for_pks('PKB76', 'PKB_probknot.ct')
#massive_test_prob_knot()
#massive_test_spot_rna()
#massive_test_ipknot()
#all_pk_programs(html)
#read_excel_overlap()
#compare_prob_to_spot()
#compare_ip_to_spot()
#check_pk_counter()
#pseudoSNitch_analysis_SPOT(html)
#nested_checker_checker()
#massive_test_spot_rna_1bp()
#massive_test_spot_rna()
def actual_pks_number():
    data_file = open('pseudobase.drct/stored_pseudobase_extractions_ALL.txt','r')
    lines = data_file.readlines()
    num_pks = 0
    for line in lines:
        for i in line:
            if i == '<' or i == '(':
                num_pks += 1
    print(num_pks)

#actual_pks_number()