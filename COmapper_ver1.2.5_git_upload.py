'''
COmapper for recombinant molecule detection from Nanopore long read sequence

Initiated on 2021.09.07.

@author: Dohwan Byun
@coworked with Jaebeom Kim

Version information 
ver0.0(210907): File input, output backbone, CIGAR transformation and genotyping
ver0.1(211118): add multiprocessing parts, recombinant read classifier
ver0.2(220422): manual input resources location, adjust read classifier with SNP number test, change output format
ver0.3(230710): change recombinant classifier setting, linux compatible
ver1.0(230725): change recombinant classifier algorithm, CO-GC detection test
ver1.1(230831): add analysis exclude regions, adjust CO-GC accuracy
ver1.2(230907): add CO frequency calculation output
ver1.2.1(230911): change CO frequency calculation, adjust CO-GC accuracy, CO site detection added, command line arguments change
ver1.2.2(231030): Filtered bases calculation correction (not include location filtered reads)
ver1.2.3(231207): command line output to output file
ver1.2.4(240412): Do not print total output
ver1.2.5(240807): git upload vesion

'''

import multiprocessing
import sys
import os
from typing import Tuple
import pandas as pd
import numpy as np
import time
import operator
from multiprocessing import Pool
import argparse
from pandas.core.frame import DataFrame
from pandas.io.parsers import read_csv, read_fwf
from functools import partial, total_ordering

class Result:
    def __init__(self, read_id, chr_num, pos, num_of_snp, snp_pos, bases, types, SNPs, length):
        self.read_id = read_id
        self.pos = pos
        self.num_of_snp = num_of_snp
        self.chr_num = chr_num
        self.snp_pos = snp_pos
        self.bases = bases
        self.types = types
        self.SNPs = SNPs
        self.length = length
        self.haplotype = self.determine_haplotype()
        self.CO_info = self.determine_CO_site()


    def determine_haplotype(self):

        string_type = self.types

        if len(string_type) < 8:
            return -1

    # filter repeat sequences
        if (self.chr_num == 1 and self.pos > 15082000 and self.pos+self.length < 15200000):
            return -1
        if (self.chr_num == 1 and self.pos > 16510000 and self.pos+self.length < 16518000):
            return -1
        if (self.chr_num == 2 and self.pos > 1000 and self.pos+self.length < 20000):
            return -1
        if (self.chr_num == 2 and self.pos > 3252000 and self.pos+self.length < 3628000):
            return -1
        if (self.chr_num == 3 and self.pos > 13588000 and self.pos+self.length < 13711000):
            return -1
        if (self.chr_num == 3 and self.pos > 14191000 and self.pos+self.length < 14224000):
            return -1
        if (self.chr_num == 4 and self.pos > 3054000 and self.pos+self.length < 3062000):
            return -1
        if (self.chr_num == 4 and self.pos > 3872000 and self.pos+self.length < 4012000):
            return -1
        if (self.chr_num == 5 and self.pos > 11184000 and self.pos+self.length < 11189000):
            return -1
        if (self.chr_num == 5 and self.pos > 11701000 and self.pos+self.length < 11739000):
            return -1
        #filter SNP wobble region
        # if (self.chr_num == 2 and self.pos > 2200000 and self.pos+self.length < 2210000):
        #     return -1
        
        if (string_type[0]==string_type[-1]):
            return 0

    # Check crossover (simple CO / COGC)
        if "CCCCLLLL" in string_type and string_type[:2]=="CC" and string_type[-2:]=="LL":
            return 1
        elif "LLLLCCCC" in string_type and string_type[:2]=="LL" and string_type[-2:]=="CC":
            return 2
        elif "LLLCLCCC" in string_type and string_type[:2]=="LL" and string_type[-2:]=="CC":
            return 11
        elif "CCCLCLLL" in string_type and string_type[:2]=="CC" and string_type[-2:]=="LL":
            return 12
        elif "LLLCCLCCC" in string_type and string_type[:2]=="LL" and string_type[-2:]=="CC":
            return 13
        elif "CCCLLCLLL" in string_type and string_type[:2]=="CC" and string_type[-2:]=="LL":
            return 14
        elif "LLLCLLCCC" in string_type and string_type[:2]=="LL" and string_type[-2:]=="CC":
            return 15
        elif "CCCLCCLLL" in string_type and string_type[:2]=="CC" and string_type[-2:]=="LL":
            return 16
        elif "LLLCCLLCCC" in string_type and string_type[:2]=="LL" and string_type[-2:]=="CC":
            return 17
        elif "CCCLLCCLLL" in string_type and string_type[:2]=="CC" and string_type[-2:]=="LL":
            return 18
        else:
            return 0


    def determine_CO_site(self): 
        string_backward = ""
        string_forward = ""
        ref_backward = ""
        ref_forward = ""

        #identify type of crossover molecule
        if self.haplotype == 1:
            ref_backward = "CCCC"
            ref_forward = "LLLL"
        elif self.haplotype == 2:
            ref_backward = "LLLL"
            ref_forward = "CCCC"
        elif self.haplotype == 11:
            ref_backward = "CLLL"
            ref_forward = "LCCC"
        elif self.haplotype == 12:
            ref_backward = "LCCC"
            ref_forward = "CLLL"
        elif self.haplotype == 13:
            ref_backward = "CCLLL"
            ref_forward = "LCCC"
        elif self.haplotype == 14:
            ref_backward = "LLCCC"
            ref_forward = "CLLL"
        elif self.haplotype == 15:
            ref_backward = "CLLL"
            ref_forward = "LLCCC"
        elif self.haplotype == 16:
            ref_backward = "LCCC"
            ref_forward = "CCLLL"
        elif self.haplotype == 17:
            ref_backward = "CCLLL"
            ref_forward = "LLCCC"
        elif self.haplotype == 18:
            ref_backward = "LLCCC"
            ref_forward = "CCLLL"
        else:
            return 0

        #left border, right border SNP identification
        for i in range(3, self.num_of_snp-3):
            t=i
            string_backward = ""
            while 1:
                if len(string_backward)==len(ref_backward) or t < 0 :
                    break
                if self.SNPs[t] == "N":
                    t-=1
                    continue
                if self.SNPs[t] != "N":
                    string_backward+=self.SNPs[t]
                    t-=1
                    continue
            if string_backward != ref_backward:
                continue

            t=i+1
            string_forward = ""
            while 1:
                if len(string_forward)==len(ref_forward) or t > self.num_of_snp-1 :
                    break
                if self.SNPs[t] == "N":
                    t+=1
                    continue
                if self.SNPs[t] != "N":
                    string_forward+=self.SNPs[t]
                    t+=1
                    continue
            if string_forward != ref_forward:
                continue

            #return CO information of read
            if string_backward == ref_backward and string_forward == ref_forward:
                lb = self.snp_pos[i]
                t=i+1
                while(1):
                    if self.SNPs[t] == "N":
                        t+=1
                        continue
                    else:
                        break
                rb = self.snp_pos[t]
                co=(lb+rb)/2
                width = rb-lb
                return [self.chr_num, lb, rb, co, width]
            else:
                continue
        return 0

    def print_result(self):
        print(self.read_id, "\t", self.chr_num, "\t", self.pos, "\t", self.num_of_snp, "\t", self.snp_pos, "\t",
              self.bases, "\t", self.types, "\t", self.SNPs, "\t", self.haplotype, "\t", self.CO_info)

    def write_result(self, df):
#        a = ([[self.read_id, self.chr_num, self.pos, self.length, self.num_of_snp, self.snp_pos, self.bases, self.types, self.Nratio, self.haplotype, self.CO_info]])
        a = ([[self.read_id, self.chr_num, self.pos, self.length, self.num_of_snp, self.snp_pos, self.bases, self.types, self.SNPs, self.haplotype, self.CO_info]])
        df = df.append(a)
        return df

def reconstruct_sequence_with_cigar(seq, cigar):

    # M: match
    # I: insertion to reference --> delete inserted bases
    # D: deletion from the reference  --> add Ns into deleted position
    # N: skipped region from the reference  --> add Ns into skipped position
    # S: soft clipping (clipped sequences present in SEQ)  --> delete from SEQ
    # H: hard clipping (clipped sequences NOT present in SEQ)  --> pass
    # P: padding (silent deletion from padded reference) --> pass
    # =: sequence match
    # X: sequence mismatch

    #threshold = 0.5

    new_seq = ""
    len_cigar = len(cigar)
    cigar_idx = 0
    seq_walker = 0
    num_non_match = 0
    len_compare = len_cigar
    while cigar_idx < len_cigar:

        # Get the number
        num = ""
        while 47 < ord(cigar[cigar_idx]) < 58:
            num = num.__add__(cigar[cigar_idx])
            cigar_idx += 1
        num = int(num)

        # Get the CIGAR code
        code = cigar[cigar_idx]
        cigar_idx += 1

        # Construct sequence
        if code == "M":
            new_seq += seq[seq_walker:seq_walker + num]
            seq_walker += num

        elif code == "I":
            seq_walker += num

        elif code == "D":
            for i in range(0, num):
                new_seq += 'N'

        elif code == "N":
            for i in range(0, num):
                new_seq += 'N'

        elif code == "S":
            seq_walker += num

            num_non_match += num
            len_compare -= num

        elif code == "H":
            num_non_match += num
            continue

        elif code == "P":
            continue

        elif code == "=":
            new_seq += seq[seq_walker:seq_walker + num]
            seq_walker += num

        elif code == "X":
            new_seq += seq[seq_walker:seq_walker + num]
            seq_walker += num

        else:
            print("something wrong")

    # Filter too many S or H (clips) containing sequence
    #if num_non_match > threshold * len_compare :
    #    new_seq = ""

    return new_seq

def CO_detect_file (pos_list, col_list, ler_list, file):
    dt = {'chr':str, 'pos':int, 'cigar':str, 'seq':str}
    f = pd.read_csv(file, delimiter='\t', names=['chr', 'pos', 'cigar', 'seq'], header=None, dtype=dt, comment='c', error_bad_lines=False)
    # Divide tsv file by chromosome number and sort

    sam_list = list()
    sam_list.append(f[f['chr'] == "Chr1"])
    sam_list.append(f[f['chr'] == "Chr2"])
    sam_list.append(f[f['chr'] == "Chr3"])
    sam_list.append(f[f['chr'] == "Chr4"])
    sam_list.append(f[f['chr'] == "Chr5"])
    for i in range(0, 5):
        sam_list[i]
        sam_list[i] = sam_list[i].sort_values(by=['pos'], axis=0)
    del f

    # Main routine #
    # Take each sequence and reconstruct according to CIGAR string,
    # and compare bases in SNP positions with the reference SNP list
    # so that determine the haplotype of each sequence

    read_cnt = 0
    col_cnt = 0
    ler_cnt = 0
    n_cnt = 0
    CO_cnt = 0
    COGC_cnt = 0
    Par_cnt = 0
    filtered_cnt = 0
    list_of_results = list()
    total_length = 0
    filtered_length = 0

    for chr_cnt in range(0, 5):  # Iterate for sequences of each chromosome
        start_idx = 0
        num_of_snp = len(pos_list[chr_cnt])
        for i, read in sam_list[chr_cnt].iterrows():  # Take each sequence
            start_pos = read['pos']
            sequence = reconstruct_sequence_with_cigar(read['seq'], read['cigar'])
            length = len(sequence)
            total_length += length
            if length == 0:
                continue
            pos_list_temp = []
            allele_temp = []
            genotype_temp = ""
            SNP_temp = ""
            snp_cnt = 0
            check_first = 0

            # Compare the bases at the SNP positions to the reference SNP list
            for j in range(start_idx, num_of_snp):
                snp_position = pos_list[chr_cnt][j]
                if 0 <= snp_position - start_pos < length:
                    if check_first == 0:
                        check_first = 1
                        start_idx = j
                    pos_list_temp.append(snp_position)
                    snp_cnt += 1
                    nucleotide = sequence[snp_position - start_pos]
                    allele_temp.append(nucleotide)
                    if nucleotide == col_list[chr_cnt][j]:
                        genotype_temp+="C"
                        SNP_temp+="C"
                        col_cnt += 1
                    elif nucleotide == ler_list[chr_cnt][j]:
                        genotype_temp+="L"
                        SNP_temp+="L"
                        ler_cnt += 1
                    else:
                        SNP_temp+="N"
                        n_cnt += 1
                elif snp_position >= start_pos + length:
                    break
            list_of_results.append(Result(i, chr_cnt + 1, start_pos, snp_cnt, pos_list_temp, allele_temp, genotype_temp, SNP_temp, length))
            read_cnt += 1

    list_of_results.sort(key=operator.attrgetter('read_id'))

    #collect CO molecule information to output dataframe
    df = pd.DataFrame()
    
    for result in list_of_results:
        if result.haplotype==1 or result.haplotype==2:
            df = result.write_result(df)
            CO_cnt += 1
            filtered_length += result.length
        elif result.haplotype > 10 :
            df = result.write_result(df)
            COGC_cnt += 1
            filtered_length += result.length
        elif result.haplotype==0:
            Par_cnt += 1
            filtered_length += result.length
        else:
            filtered_cnt += 1

    #calculate CO frequency for each files
    
    print("File name: ", file)
    print("Number of sequences: ", read_cnt)
    print("Total length: ", total_length)
    print("Filtered length: ", filtered_length)
    print('Col SNP count: ',col_cnt)
    print('Ler SNP count: ', ler_cnt)
    print('Not determined SNP count:', n_cnt)
    print('Total SNP number: ', (col_cnt+ler_cnt+n_cnt))
    print("CO molecule: ", CO_cnt)
    print("COGC molecule: ", COGC_cnt)
    print("Parental molecule: ", Par_cnt)
    print('filtered molecule: ',filtered_cnt)
    print('Total CO: ', (CO_cnt+COGC_cnt))
    if col_cnt+ler_cnt+n_cnt != 0:
        snp_density = total_length/(col_cnt+ler_cnt+n_cnt)
        print('Average SNP density: ', snp_density)
        if read_cnt != 0:
            print("Average SNPs per read: ", (col_cnt+ler_cnt+n_cnt)/read_cnt)
    if min(col_cnt, ler_cnt) != 0:
        cl_ratio = (col_cnt+ler_cnt)/min(col_cnt, ler_cnt)/2
        print("C/L ratio calibration: ", cl_ratio)
    if filtered_length != 0:
        CO_per_Gb = (CO_cnt+COGC_cnt)/filtered_length*1000000000*cl_ratio
        print('Crossovers per Gb: ', CO_per_Gb)
    if Par_cnt != 0:
        CO_per_Par = (CO_cnt+COGC_cnt)/Par_cnt*100*cl_ratio
        print("Crossovers per Parental: ", CO_per_Par)

    return(df)

def main():
    time0 = time.time()
    INPUTFOLDER = ""
    SNPLIST = ""
    CPUNUMBER = 0
    OUTPUTFILE = ""

    try:
        #command line arguments
        parser = argparse.ArgumentParser(description='process input folder location')
        parser.add_argument('--inputfolder', '-i', help = 'input folder location, default: current location', default = os.getcwd())
        parser.add_argument('--snpfile', '-s', help='snp file location, default: collerF2.masked.tiger.txt', default=os.getcwd()+'/collerF2.masked.tiger.txt')
        parser.add_argument('--threads', '-t', help = 'thread number, default: 6', default = 6)
        parser.add_argument('--output', '-o', help = 'output file name, default: result.csv', default = "result.csv")
        args = parser.parse_args()
        INPUTFOLDER = args.inputfolder
        SNPLIST = args.snpfile
        CPUNUMBER = int(args.threads)
        OUTPUTFILE = args.output
    except:
        print("input error -h for help")
        sys.exit(2)

    # Read a file of reference SNPs
    #manual input folder and snp list
    # INPUTFOLDER = '/datasets/data_2/dohwan/Nanopore/231123_randomsampling/recq4_pollen/output'
    # SNPLIST = '/datasets/data_1/dohwan/Nanopore/resources/collerF2.masked.tiger.txt'
    # OUTPUTFILE = "result_recq4_repeat_ver1_2_3.csv"
    # CPUNUMBER=10
    
    print("input folder:", INPUTFOLDER)
    print("snpfile: ", SNPLIST)
    print("threads: ", CPUNUMBER)
    print("output file name: ", OUTPUTFILE)
    snp_all = pd.read_csv(SNPLIST, delimiter='\t')

    # Divide reference SNPs by chromosome number
    snp_list = list()
    for i in range(1, 6):
        snp_list.append(snp_all[snp_all['Chr'] == i])
    del snp_all

    # Make lists of (SNP position and Col/Ler) for each chromosome
    pos_list = list()
    col_list = list()
    ler_list = list()
    for i in range(0, 5):
        pos_list.append(np.array(snp_list[i]['Pos']))
        col_list.append(np.array(snp_list[i]['Col']))
        ler_list.append(np.array(snp_list[i]['Ler']))
    del snp_list

    file_list = os.listdir(INPUTFOLDER)

    #manual file location
    file_list_t = []
    for file in file_list :
        file_path = os.path.join(INPUTFOLDER, file)
        file_list_t.append(file_path)
    print(file_list_t)

    # Read tsv input file, detect CO
    df = pd.DataFrame()
    pool = multiprocessing.Pool(processes=CPUNUMBER)
    func = partial(CO_detect_file, pos_list, col_list, ler_list)
    df = df.append(pool.map(func, file_list_t), ignore_index=True)
    pool.close()
    pool.join()

    # create output file
    if(df.empty==False):
        df.columns=(['read_id', 'chr_num', 'pos', 'length', 'num_of_snp','snp_pos', 'bases', 'types', 'SNP_genotype', 'haplotype', 'CO_info'])
        df.to_csv(OUTPUTFILE, mode='w')

    else:
        print("No CO detected")
    #print('total reads: ',total_read)
    print(time.time() - time0)
    
if __name__ == '__main__':
    main()
