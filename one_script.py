#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#sliding window PLEASE READ!!!!
# arguments: -f for filename; -from -to for range; -ppm if to provide accuracy in ppm; -Da for acuracy in daltons; -w for window_size; -m for mode of search; -p for pattern of interest
# -m and -p are not required if there is no pattern of interest.

import argparse, re

value_list = []
bins = {}
aa_score_list = []
error_list = ['B', 'J', 'O', 'U', 'X', 'Z']

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file',
    help='chose the file', required=True)
parser.add_argument('-from', '--from_range',
    help='chose the range. From:', default=1000)
parser.add_argument('-to', '--to_range',
    help='chose the range. To:', default=1500)
parser.add_argument('-ppm', '--parts_per_million',
    help='enter the value for ppm')
parser.add_argument('-Da', '--Da_accuracy',
    help='enter the value for mass accuracy in Da', default=0.2)
parser.add_argument('-w', '--window_size',
    help='define the window size', required=True)
parser.add_argument('-m', '--mode',
    help='chose the mode of search [any, start, end]')
parser.add_argument('-p', '--pattern',
    help='define the pattern ')


args = parser.parse_args()
window_size = float(args.window_size)
x = float(args.from_range)
y = float(args.to_range)
Da = args.Da_accuracy

def search_any(input_data, pattern):
    print ('Peptide ions containing {0}:'.format(pattern))
    for line in input_data:
        line = line.rstrip("\n")
        f = line.split()
        seq = f[5]
        m_aa = re.search(r'{0}'.format(pattern), seq)
        if m_aa:
            aa_score_list.append(line)

            
def search_start(input_data, pattern):
    print ('Peptide ions starting with {0}:'.format(pattern))
    for line in input_data:
        line = line.rstrip("\n")
        f = line.split()
        seq = f[5]
        m_aa = re.search(r'^{0}'.format(pattern), seq)
        if m_aa:
            aa_score_list.append(line)


def search_end(input_data, pattern):
    print ('Peptide ions ending with {0}:'.format(pattern))
    for line in input_data:
        line = line.rstrip("\n")
        f = line.split()
        seq = f[5]
        m_aa = re.search(r'{0}$'.format(pattern), seq)
        if m_aa:
            aa_score_list.append(line)

try:
    data = open(args.file).readlines()    #Error catch
except IOError:
    print ('Could not open file',args.file)
    exit()

if float(args.from_range) == 0 and args.parts_per_million:  #Error catch
    print ('You cannot use range from 0 with ppm! Try using a number close to null instead')
    exit()

if args.mode and args.pattern:
    aa = args.pattern.upper()
    for element in error_list:
        match_res = re.search(r'{0}'.format(element), aa)
        match_res1 = re.search(r'[0-9]+',aa)
    if match_res:
         print('Please enter the correct residues. Amino acids are not represented by {0}'.format(error_list))
         exit()
    elif match_res1:
         print('Residues cannot be represented by numbers!')
         exit()

if args.mode == 'any':
    search_any(data, aa)
elif args.mode == 'start':
    search_start(data, aa)
elif args.mode == 'end':
    search_end(data, aa)
            

def range_(data, rfrom, rto, list_):
    for i in data:
        f = i.split()
        z = float(f[2])
        if z >= rfrom and z <= rto:
            list_.append(z)

if args.mode and args.pattern and aa_score_list == []:
    print ('NO MATCHES WITH THE PATTERN')
    exit()
elif args.mode and args.pattern:
    range_(aa_score_list, x, y, value_list)
else:
    range_(data, x, y, value_list)

if value_list == []:
    print ('No values in the given mass range!')  #Error catch
    exit()
#print (value_list, len(value_list)) #checkpoint
        
s = sorted(value_list)

m2 = x + window_size
m1 = m2 - window_size
while m2 <= y:
    for i in s:
        if m2 not in bins.keys():
            bins[m2] = 0
        if i >= m1 and i <= m2:
            bins[m2] += 1
    if args.parts_per_million:
        step = float(args.parts_per_million) * m1 / 10**6
        m1 += step
        m2 += step
    elif Da and float(Da) <= window_size:
        m1 += float(args.Da_accuracy)
        m2 += float(args.Da_accuracy)
    elif float(Da) > window_size:          #Error catch
        print ('The mass accuracy (step) is larger than the window size!')
        exit()

sorted_k = sorted(bins.keys())
for i in sorted_k:
    print (i, bins[i])

print ('With range from {0} to {1}, and selected window size {2}'.format(x, y, window_size))
print ('There are {0} values in this range of mass/charge:'.format(len(s)))
if args.parts_per_million:
    print ('The ppm entered was {0}'.format(args.parts_per_million))
elif Da:
    print ('The step entered in daltons was: {0}'.format(Da))
