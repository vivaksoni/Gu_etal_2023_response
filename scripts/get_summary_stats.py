#script takes .ms file and calculates summary statistics across sliding windows using libsequence

#python3  stats_sliding_window.py -msFile "rep1.ms" \
#-fixedFile "rep1.fixed" \
#-outPath  "rep1.stats" \
#-regionLen 30000 -samples 100 -l_threshold 0.025 -u_threshold 0.5

from __future__ import print_function
import libsequence
import sys
import pandas as pd
import math
import numpy as np
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-msFile', dest = 'msFile', action='store', nargs = 1, type = str, help = 'path to ms file (.ms format)')
parser.add_argument('-outPath', dest = 'outPath', action='store', nargs = 1, type = str, help = 'path to output files')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of simulated region')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples in .ms file')
parser.add_argument('-l_threshold', dest = 'l_threshold', action='store', nargs = 1, type = float, help = 'AF minimum threshold')
parser.add_argument('-u_threshold', dest = 'u_threshold', action='store', nargs = 1, type = float, help = 'AF maximum threshold')

args = parser.parse_args()
chr_len =  args.regionLen[0]
samples = args.samples[0]
l_threshold = args.l_threshold[0]
u_threshold = args.u_threshold[0]
ms_file = args.msFile[0]
out_path = args.outPath[0]

#function reads through .fixed file, 
def read_fixed_mutations(f_fixed, chr_len, N, burnIn):
    #Create empty dict to store substitutions
    d_subs = {}
    #Loop through lined in.fixed file
    for line in f_fixed:
        #strip newline character
        line1 = line.strip('\n')
        #split line into list
        lst = line1.split()
        #Skip first two lines to move straight to data lines
        if line1[0]!="#" and lst[0]!="Mutations:" and int(lst[-1])>(N*burnIn):
            #Estimate position using decimal notation
            posn = float(lst[3])/float(chr_len)
            #Assign base position as key
            #With each line, if position occurs again, value is incremented
            #(To account for repeat mutation in the same position)
            d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs 



def avg_divergence_win(d_subs, start, end):
    s_sum = 0
    for posn in d_subs.keys():
        if float(posn) <= end and float(posn) > start:
            s_sum = s_sum + 1
    return s_sum



#Get number of segregating sites from .ms file
def get_S(f_ms):
    #Return no. of segregating sites (line 1 of .ms file)
    pos_lines = [1]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            S = line.split()[1]
    #Set file pointer to start of file
    f_ms.seek(0)
    return S



#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, samples):
    l_Pos = [] #list of positions of SNPs
    l_Genos = [] #list of alleles
    d_tmp = {} #dict to store individual allele info for each individual (values) at each site (keys)

    #positions on line 2
    pos_lines = [2]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            #store positions in list
            pos_list  = line.split()
    
    #Set file pointer to start of file
    f_ms.seek(0)
    i = 0
    #Loop through positions, storing in list
    for pos in pos_list[1:]:
        #Append position to l_Pos (after converting to float)
        l_Pos.append(float(pos))
        #Add dictionary key for each position, with empty value
        d_tmp[str(i)] = ""
        i += 1  

    #genotypes on line 3 onwards (use samples argument to determine length of file)
    g_lines = [x for x in range(3, samples + 4)]
    #Loop through lines (ie individuals)
    for position, line in enumerate(f_ms):
        if position in g_lines:
            #Remove newline character
            line1 = line.strip('\n')
            i = 0
            #For each individual, loop through each site, appending allele information for that individual to 
            #the site number in dict
            while i < len(line1):
                d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                i = i + 1

    f_ms.seek(0)

    #Create nested list of positions and genotypes
    l_data = [[j, d_tmp[str(i)]] for i,j in enumerate(l_Pos)]
    return(l_data)


#Function to filter SNPs by a minimum allele frequency
def filter_by_MAF(l_data, l_threshold, u_threshold, samples):
    res = []
    masked = []
    for i in l_data:
        if(((i[1].count('1')/samples)>l_threshold) & ((i[1].count('1')/samples)<u_threshold)):
            res.append(i)
        else:
            masked.append(i)
    return([res, masked])


#Function to create dictionary of polySIM summary statistics. Returns dictionary of summary stats
def get_polySIM_stats(sd):
    #Create polysim object
    ps = libsequence.PolySIM(sd)
    #Create list of methods (ie polySIM summaryStats)
    a = [method for method in dir(ps) if callable(getattr(ps, method)) if not method.startswith('_')]
    #Loop through methods, storing names as keys, and estimates as values in dict
    ss_dict = {}
    for method in a:
        ss_dict[method] = getattr(ps, method)()
        
    return(ss_dict)



#Function to create dictionary of LD stats. Returns dictionary.
def get_LD_stats(sd):
    ld = libsequence.ld(sd)
    df = pd.DataFrame(ld)
    ss_dict = {}
    ss_dict['meanrsq'] = sum(df['rsq'])/len(df['rsq'])
    ss_dict['meanD'] = sum(df['D'])/len(df['D'])
    ss_dict['meanDprime'] = sum(df['Dprime'])/len(df['Dprime'])
    return(ss_dict)


#Read in .ms file, create sd object for libsequence
f_ms = open(ms_file, 'r')
S = get_S(f_ms)
l_data = get_nested_data_list(f_ms, 100)
#Filter data via allele frequency threshold
filtered = filter_by_MAF(l_data, l_threshold, u_threshold, samples)
l_data_filtered = filtered[0]
l_data_masked = filtered[1]

#Convert to df, and convert positions into integer coordinates
fdf = pd.DataFrame(l_data_filtered, columns=['position', 'genotype'])
fdf['pos'] = [int(x)+1 for x in (np.round(fdf['position']*chr_len, 0))]
#Determine if site is non-synonymous or synonymous
fdf['csq'] = np.where(fdf.pos%3==0, 'S', 'NS')
#Create nested lists for synonymous and non-synonymous sites
S = fdf[fdf.csq=='S'][['pos', 'genotype']].values
NS = fdf[fdf.csq=='NS'][['pos', 'genotype']].values

#Repeat for masked data
fdf2 = pd.DataFrame(l_data_masked, columns=['position', 'genotype'])
fdf2['pos'] = [int(x)+1 for x in (np.round(fdf2['position']*chr_len, 0))]
#Determine if site is non-synonymous or synonymous
fdf2['csq'] = np.where(fdf2.pos%3==0, 'S', 'NS')
#Create nested lists for synonymous and non-synonymous sites
S2 = fdf2[fdf2.csq=='S'][['pos', 'genotype']].values
NS2 = fdf2[fdf2.csq=='NS'][['pos', 'genotype']].values

sd_S = libsequence.SimData(S)
sd_NS = libsequence.SimData(NS)

d_NS = get_polySIM_stats(sd_NS)
d_S = get_polySIM_stats(sd_S)
#Convert dict to dataframe and append to file
df = pd.DataFrame.from_dict(d_NS, orient='index').T
df['num_sites'] = len(NS2)
df.to_csv(out_path + '_NS.txt', sep='\t', index=False, header=False, mode='a')
df = pd.DataFrame.from_dict(d_S, orient='index').T
df['num_sites'] = len(S2)
df.to_csv(out_path + '_S.txt', sep='\t', index=False, header=False, mode='a')
print ("Stats output to " + out_path)




