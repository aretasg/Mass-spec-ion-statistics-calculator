#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#PLEASE READ!!!!
# arguments: -f for filename; -from -to for range; -ppm to provide accuracy in ppm; -Da for acuracy in daltons; -w for window_size; -m for mode of search; -p for pattern of interest
# -m and -p are not required if there is no pattern of interest.

import argparse     #importing modules
import re
import pylab

value_list = [] #a list of values in the specified range
bins = {} #a dictionary for counting accurances of m/z values in the specified window size
seq_score_list = []  #a list of sequences that have pattern of interest 
error_list = ['B', 'J', 'O', 'U', 'X', 'Z']  #letters not used to mark amino acids

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file',
    help='chose the file', required = True)  #a flag to specify the input file
parser.add_argument('-from', '--from_range',
    help='chose the range. From:', default = 1000) #a flag to specify the starting point of the range; default is 1000 m/z
parser.add_argument('-to', '--to_range',
    help='chose the range. To:', default = 1500)    #a flag to specify the ending point of the range; default is 1500 m/z
parser.add_argument('-ppm', '--parts_per_million',
    help='enter the value for ppm')                 #a flag to specify accuracy as parts per a million
parser.add_argument('-Da', '--Da_accuracy',
    help='enter the value for mass accuracy in Da', default = 0.2)  #a flag to specify accuracy in Daltons; default is 0.2 Da
parser.add_argument('-w', '--window_size',
    help='define the window size', required = True)  #a flag for window size for the sliding window
parser.add_argument('-m', '--mode',
    help='chose the mode of search [any, start, end]') #a flag for the type of pattern search to be performed
parser.add_argument('-p', '--pattern',
    help='define the pattern')  #a flag for the patern of interest to search e.g. CAGIY (one letter code for amino acid)


args = parser.parse_args()
window_size = float(args.window_size)  #a variable for the window size as provided from CL
start_point_range = float(args.from_range)  #a variable for the starting point of the range as provided from CL
end_point_range = float(args.to_range)    #a variable for the ending point of the range as provided from CL
Da = args.Da_accuracy       #a variable for the accuracy in Daltons as provided from CL
input_file = args.file      #a variable for the name of the input file

#functions###
def search_any(input_data, pattern): #a function to search data for patterns. Starting with/ending with or just containing them anywhere
    print ('Counting peptide ions containing the residue pattern: {0}'.format(pattern)) #print message informing the specified pattern from CL
    for line in input_data:             #iterates the input file
        line1 = line.rstrip("\n")        #removes the next line sign
        f = line1.split()            #splits the line
        seq = f[5]              #save the fifth element/sequence of the line as variable seq
        m_aa = re.search(r'{0}'.format(pattern), seq) #search for the pattern in the sequence and save it to variable m_aa
        if m_aa:                        #if there is a hit 
            seq_score_list.append(line1)  #add the sequence to the hit list

            
def search_start(input_data, pattern): #a function to search data for sequences starting with the residues provided with -p
    print ('Counting peptide ions starting with {0}:'.format(pattern))
    for line in input_data:
        line1 = line.rstrip("\n")
        f = line1.split()
        seq = f[5]
        m_aa = re.search(r'^{0}'.format(pattern), seq)  #looks for the sequences starting with the given pattern
        if m_aa:
            seq_score_list.append(line1)


def search_end(input_data, pattern): #a funtion to search data for sequences ending with the residues provided with -p
    print ('Counting peptide ions ending with {0}:'.format(pattern))
    for line in input_data:
        line1 = line.rstrip("\n")
        f = line1.split()
        seq = f[5]
        m_aa = re.search(r'{0}$'.format(pattern), line)  #looks for the sequences ending with the given pattern
        if m_aa:
            seq_score_list.append(line1)
            
def range_function(sequences, start_point, end_point, range_list):  #a function to iterate over the input; compare m/z values to a specified range (-form/-to)
    for i in sequences:                    #iterate over the list containing sequences
        f = i.split()                       #split the line
        seq = float(f[2])                     #save third element/ m/z value 
        if seq >= start_point and seq <= end_point: #checks if the value is in the range
            range_list.append(seq)      #adds the value to the list for m/z values in the specified range

#################end of functions
try:                                        #attempts to open and read the input
    data = open(input_file).readlines()
except IOError:
    print ('Could not open file', input_file)  #Displays an error message if unsuccessfull
    exit()

if start_point_range < 0 or end_point_range < 0:  #check if the range value is not negative; m/z values cannot be negative
    print ('m/z range values cannot be negative')
    exit()
    
if float(start_point_range) == 0 and args.parts_per_million:  #Checks if range starts with 0 when ppm accuracy is selected
    print ('You cannot use range starting with 0 together with the ppm option! Try using a number close to null instead')
    exit()

if args.mode and args.pattern:        #Checks if the pattern is sensible and can be used to search the data
    pattern_upper_case = args.pattern.upper() #makes pattern uppercase; so the user can type lower case in the CL
    for element in error_list:
        match_res = re.search(r'{0}'.format(element), pattern_upper_case)  #checks if the pattern has any letters not representing amino acids
        match_res1 = re.search(r'[0-9]+', pattern_upper_case)  #checks if the pattern has any numbers (amino acids cannot be represented as numbers
    if match_res:                           #if the pattern has any letters in the error_list; displays the message 
         print('Please enter the correct residues. Amino acids are not represented by {0}'.format(error_list))
         exit()
    elif match_res1:  #displays the message if the patter contains any numbers
         print('Residues cannot be represented by numbers!')
         exit()

        
#Depending on the preference specified in CLI, one out of three pattern search options are selected and
#a respective function is called
if args.mode == 'any': #checking the mode of search
    search_any(data, pattern_upper_case) #running a function search_any
elif args.mode == 'start':
    search_start(data, pattern_upper_case)
elif args.mode == 'end':
    search_end(data, pattern_upper_case)

    
if args.mode and args.pattern and seq_score_list == []:  #prints a message if no sequences match the pattern of interest
    print ('NO MATCHES WITH THE PATTERN: ' + string(pattern_upper_case))    
    exit()              #and closes the program
elif args.mode and args.pattern:  #select values in the pattern hit list that are in range
    range_function(seq_score_list, start_point_range, end_point_range, value_list)   #adds m/z values of pattern hits that are in the range
else:
    range_function(data, start_point_range, end_point_range, value_list)  #adds m/z values in range to the value_list if no pattern was specified


if value_list == []:        ##Prints a warning message if there are no values in the m/z range specified and exits the program
    print ('No values in the given mass range!')
    exit()
#print (value_list, len(value_list)) #checkpoint

####Sliding window
sorted_list = sorted(value_list)  #sorts the value_list

window_end = start_point_range + window_size    #the end point of the window
window_start = window_end - window_size         #the start point of the window

while window_end <= end_point_range:           #while the end of window is less or equals to the end point of the range
    for i in sorted_list:                      #for every m/z value in the sorted list (value_list)
        if window_end not in bins.keys():       #creates a key for the end position of the window in the bins dictionary
            bins[window_end] = 0                #set the value of the key as zero
        if i >= window_start and i <= window_end: #if the m/z value is in the position of the window; count the occurence
            bins[window_end] += 1
    if args.parts_per_million:                              #if ppm was selected the step to move the window is counted appropriately
        step = float(args.parts_per_million) * window_start / 10**6  #the step with the accuracy applied
        window_start += step
        window_end += step
    elif Da and float(Da) <= window_size:       #move the window based on accuracy provided in Daltons
        window_start += float(args.Da_accuracy)
        window_end += float(args.Da_accuracy)
    elif float(Da) > window_size:          #check if the window size is larger than the accuracy in Daltons
        print ('The mass accuracy (step) is larger than the window size!')
        exit()

sorted_keys = sorted(bins.keys()) #sort the bins dictionary by the key; from lowest to highest
for i in sorted_keys:   #prints the results the the end point of window positions followed by the number of m/z values in that position
    print (i, bins[i])

 ##prints some information about the input file; the pattern; accuracy selected
print ('In the range from {0} to {1} m/z and window size {2};'.format(start_point_range, end_point_range, window_size))
print ('There are {0} values in this range of mass/charge'.format(len(sorted_list)))
if args.parts_per_million:
    print ('The parts per a million (ppm) accuracy entered was {0}'.format(args.parts_per_million))
elif Da:
    print ('The accuracy entered in Daltons was: {0}'.format(Da))
if args.mode and args.pattern:
    print ('The pattern entered was: {0}'.format(pattern_upper_case))
    
#creates a plot and saves it in the working directory to allow graphical visualisation of the results

pylab.figure(1)
pylab.bar(list(sorted(bins.keys())), bins.values(), width=0.1, linewidth=2,
        edgecolor='black')
pylab.ylabel ("Intensity, rel. units")
pylab.xlabel ("m/z, Th")

if args.mode and args.pattern:
    pylab.title('MS_spectra for {0}; Sequence: {1}; {2}-{3} m/z'.format(input_file, pattern_upper_case, start_point_range, end_point_range))
    pylab.savefig('MS_spectra_{0}_SEQ_{1}.png'.format(input_file, pattern_upper_case))
else:
    pylab.title('MS_spectra for {0}; {1}-{2} m/z'.format(input_file, start_point_range, end_point_range))
    pylab.savefig('MS_spectra_{0}.png'.format(input_file))
#plt.show()

