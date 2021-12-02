#SNPfold_ACI.py
#this script is for ACI usage and serves as a work-around so that the .pbs file can switch python versions while running
#many of functions here are derived from blast_output_reader.ACI
import subprocess
import os
import time
from  blast_output_reader_ACI import SNPfold_significance_multiple_SNPs as sfsms
from  blast_output_reader_ACI import SNPfold_significance as sfs
from blast_output_reader_ACI import store_real_riboSNitches as rsrSN
import argparse


home_directory_path = '/storage/work/p/pcf5065/'

parser = argparse.ArgumentParser()

parser.add_argument('-i', type=str, help='input file directory')

args = parser.parse_args()

input_dir = (args.i)
os.chdir(input_dir.replace('/',''))

#this is the intermediate file generated in the primary RSC module script
file = open('intermediate.txt', 'r')

#finds the correct file from the linked RSC part 1 script that contains SNP data
def find_all_snps(this_run_num):
    right_file_name = ''
    for file in os.listdir(os.getcwd()):
        if ('all_snps_' in file) and (file.split('all_snps_')[1].split('_')[0].isdigit() == True):
            #finds the file with the same run ID
            if int(this_run_num) == int(file.split('all_snps_')[1].split('_')[0]):
                right_file_name = file
                break
    if right_file_name == '':
        for file in os.listdir('collected_snps_Pass'):
            if ('all_snps_' in file) and (file.split('all_snps_')[1].split('_')[0].isdigit() == True):
                if int(this_run_num) == int(file.split('all_snps_')[1].split('_')[0]):
                    right_file_name = 'collected_snps_Pass/'+file
                    break
    #returns the file with the SNP data
    return right_file_name

#adds SNPfold data to master SNP data set
def add_results_to_all_snps(snp_list, lines, snp_string, position, results):
    if 'Os' in position and '_' not in position and '-' not in position:
        #print(position)
        #position = position[15:]
        #pos_num = position[1]
        #for i in position:
        #    if i.isdigit() == True:
        #        pos_num += i
        #position = int(pos_num)
        position = int(snp_list[0])
    elif '_' in position or '-' in position:
        return None
    else:
        #pos_num = snp_string[1]
        #for i in snp_string:
        #    if i.isdigit() == True:
        #        pos_num += i
        position = int(snp_list[2])
    wt_base = str(snp_string[0])
    mu_base = str(snp_string[-1])
    output_string = ''
    for num_line, line in enumerate(lines):
        output_line = line
        data = line.split('\t')
        if len(data) > 2:
            if (str(data[0]) == wt_base) and (int(data[1]) == position) and (str(data[2]) == mu_base):
                if data[7].isdigit() == True:
                    rna_snp=int(data[7])
                    if results[0] == True:
                        rna_snp += 1
                        output_line = '\t'.join(data[:7])+'\t'+str(rna_snp)+'\t'+'\t'.join(data[8:])
                        #output_line = line.split('rbSN: ')[0]+str(rna_snp)+line.split('rbSN: ')[1][1:]
        output_string += output_line
    return output_string

#rewrites a fasta file without the head for SNPfold input
def clean_fasta(input_file_name_clean):
    input_file_clean = open(input_file_name_clean, 'r')
    input_lines_clean = input_file_clean.readlines()
    new_seq_clean = ''
    if input_lines_clean[0][0] == '>':
        for line in input_lines_clean[1:]:
            new_seq_clean += line
        output_file_name_clean = input_file_name_clean.replace('.fasta', '_clean.fasta')
        output_file = open(output_file_name_clean, 'w')
        output_file.write(new_seq_clean)
        output_file.close()
        return output_file_name_clean
    else:
        return input_file_name_clean

def run_snpfold(snp_list, snp_string_input, snp_5_prime, pk_filename, input_file, is_RNAsnp_sig, dir_name, this_run_num, default='True', output_dir='../../putative_SNitches/'):
    #formatting for file naming clarity
    output_dir = home_directory_path+output_dir
    if default == 'True':
        input_ID = input_file.split('/')[1][:-11]
    else:
        input_ID = input_file
    #configuring the SNP from RNAsnp to SNPfold
    snp_string = snp_string_input.replace('_', ',')
    snp_string = snp_string_input.replace('-', ',')
    insert = ''
    #This ensures that the snp is singular, not two SNPs at once. If it is singular, a p value is calucated
    #if ',' not in snp_string:
    #    insert = ' -A'
    #the following block of code creates output files and runs SNPfold
    SNPfold_output_name = 'SNPfold_output%s_%s.txt' % (snp_string, input_ID)
    SNPfold_output = open(SNPfold_output_name, 'w')


    os.chdir('../programs/SNPfold-1.01/')
    #modifies the sequence file to conform to SNPfold
    snpfold_input_file_name = clean_fasta(output_dir+pk_filename)
    #runs SNPfold
    python_2_command = 'python2.7 SNPfold_commandline.py%s %s %s > %s' % (insert, snpfold_input_file_name, snp_string, output_dir+SNPfold_output_name)
    os.system(python_2_command)
    os.chdir(output_dir)
    SNPfold_output = open(SNPfold_output_name, 'r')
    #this counters an error that was being encountered due to a corrupted file
    test_lines = SNPfold_output.readlines()
    new_lines = ''
    for test in test_lines:
        if 'AA' not in test:
            new_lines += test
    #the next block of code counters the corrupted files
    new_filename = 'SNPfold_outputfile%s_%s.txt'%(snp_string, input_ID)
    new_file = open(new_filename, 'w')
    new_file.write(new_lines)
    new_file.close()
    #if ',' in snp_string:
    #uses imported function to read SNPfold output data
    is_SNPfold_sig = sfsms(SNPfold_output_name, is_RNAsnp_sig)
    #stores true riboSNitches in separate file
    rsrSN(snp_5_prime, snp_string, is_SNPfold_sig, new_lines, input_file, False)
    try:
        os.remove(new_filename)
    except:
        pass
    #this block of code edits the SNP output file
    snp_string = snp_string.replace(',','-')
    #locates the SNP filename
    all_snps_file_name = find_all_snps(this_run_num)
    all_snps_file = open(all_snps_file_name, 'r')
    all_snps_lines = all_snps_file.readlines()
    #adds the new data regarding the SNPs to the master data file
    new_lines = add_results_to_all_snps(snp_list, all_snps_lines, snp_string, snp_5_prime, is_SNPfold_sig)
    all_snps_file.close()
    all_snps_file = open(all_snps_file_name, 'w')
    all_snps_file.write(new_lines)

#function that reformats a string-coerced list back into a real list
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

#this is the beginning of the SNPfold portion of the RSC module
print('\n\nBEGIN SNPFOLD')

lines = file.readlines()
#loops through all the SNPs
for line in lines:
    if line == ' ':
        break
    #collects the SNP data
    data_items = line.split('\t')
    snp_string_input = data_items[1]
    snp_5_prime = data_items[0]
    pk_filename = data_items[2]
    input_file = data_items[3]
    is_RNAsnp_sig = data_items[4]
    this_run_num = int(data_items[5])
    #converts the string into a real list
    snp_list = construct_snp_list(data_items[6])
    dir_name = data_items[7]
    default = data_items[8].replace('\n', '')
    #print(default, 'iiiii')
    #principle function of the script
    run_snpfold(snp_list, snp_string_input, snp_5_prime, pk_filename, input_file, is_RNAsnp_sig, dir_name,this_run_num,default, input_dir)
#this block removes and moves files
os.remove('intermediate.txt')    
for file in os.listdir(os.getcwd()):
    if 'ref_' in file and 'NC' not in file and 'putative' not in file:
        os.remove(file)
    if 'SNPfold' in file:
        try:
            os.rename(file, 'putative_SNitches/riboSNitch_programs_data/'+file)
        except:
            os.rename(file, 'riboSNitch_programs_data/'+file)
#for file in os.listdir(os.getcwd()):
#    if 'SNPfold' in file or 'Os' in file or 'RNAsnp' in file:
#        os.remove(file)
print(time.ctime(time.time()))