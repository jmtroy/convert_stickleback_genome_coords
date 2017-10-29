#!/usr/bin/python

# convert_gasAcu1S1_to_gasAcu1.py
# Modifications
# 05-28-2017 Joe Troy Created

# A linux/MAC OS command line program to convert stickleback fish GasAcu1 genome
# coordinates to GasAcu1S1 genome coordinates.  Assumes a tab delimited file and writes
# output to standard out

# for example: convert chromosome name from groupI to chrI, groupII to chrII, etc..
# and also convert scaffold_37 1000  2000   to chrUn 5107261 5106261 
# (per the ChrUn_to_scaffolds.txt conversion table)

# Outline of this code

# 1) accept arguments for 
#	a) file to convert 
#	b) column with chromosome 
#	c) column with start
#	d) column with end
#	e) name of conversion file, default is ChrUn_to_scaffolds.txt is current working directory
#	f) column of Scaffold names, default is column 1
#	g) column of chrom name to convert to, default is column 3
#	h) column of start position in new (converted to) chrom name, default is column 4
#	i) number of header rows in the input file (a), that will not be converted, default is 0

# 2) read in conversion file, and create a pandas data frame with the data

# 3) read in file to convert 
#		3a) loop through each row and convert chromosome, start and end to new values
#			conversion from groupI to chrI, groupII to chrII, etc.. is done with a simple
#			replace and the start and end are not changed in this case
#			scaffold_xxxx to chrUn conversion are done using the conversion table
#			and in this case the start and ends are updated.
#
#			After conversion the new rows (including the header rows) are sent to standard out (stdout)

# Note: below are the steps I used to install the pandas python library on my Mac Book
# $sudo easy_install pip
# $sudo pip install pandas

# import libraries

import glob
import os
import sys
import argparse
import time
import subprocess
import pandas as pd

# setup argparse (set up command line arguments and help)

parser = argparse.ArgumentParser(description='''convert_gasAcu1S1_to_gasAcu1.py is a linux/MAC OS command
line program to convert files with stickleback gasAcu1S1 gemome coordinates to
gasAcu1 genome coordinates. for example: convert chromosome name from groupI to chrI, groupII to chrII, etc..
And also convert scaffolds to chrUn. for example: scaffold_37 1000  2000   to chrUn 5107261 5108261.

To know how to convert scaffolds to chrUn, the program requires a file called ChrUn_to_scaffolds.txt
								
run the program by at the command line: 
python convert_gasAcu1S1_to_gasAcu1.py --input_file {tab file to convert} --chrom_col {column #} --start_col {column #} 
--end_col {column #} --header_rows {number of header rows in input file}. 
This runs the program with the defaults, see the other options for more information. ''')

parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--input_file', nargs='?' , dest = 'input_file', required=True, help='Specify the tab separated input file to convert (output is sent to stdout --REQUIRED--')
parser.add_argument('--chrom_col', nargs='?' , dest = 'chrom_col', type=int, required=True, help='column with chromosome name --REQUIRED--')
parser.add_argument('--start_col', nargs='?' , dest = 'start_col', type=int, required=True, help='column with start position --REQUIRED--')
parser.add_argument('--end_col', nargs='?' , dest = 'end_col', type=int, required=True, help='column with end position --REQUIRED--')
parser.add_argument('--header_rows', nargs='?' , dest = 'header_rows', type=int, default = 0, help='number of header rows, default is 0')
parser.add_argument('--conversion_file', nargs='?' , dest = 'conversion_file', default = 'ChrUn_to_scaffolds.txt', help='Specify the tab separated conversion (mapping) file.')
parser.add_argument('--cf_frm_name_col', nargs='?' , dest = 'cf_frm_name_col',  default = 'Scaffold', help='Within the conversion file, the column name with scaffold name')
parser.add_argument('--cf_to_name_col', nargs='?' , dest = 'cf_to_name_col',  default = 'Chr', help='Within the conversion file, the column name with new chrom name, chrUn for example')
parser.add_argument('--cf_start_col', nargs='?' , dest = 'cf_start_col',  default = 'Start', help='Within the conversion file, the column name with new start position to be used as the offset in the conversion')

args = parser.parse_args()

# read in conversion file using panda to create a dataframe
convert_df = pd.read_csv(args.conversion_file,sep='\t')
# print(convert_df.head(3))  ## this was used as a test

# below we read in each line of the input file
#	send the line to standard out "as is" if it is a header row
#	convert group... to chr... if the chromosome starts with 'group', send to standard out
#	if does not begin with group, use the convert_df to find new name and update start and end, send to standard out
row_counter = 0
with open(args.input_file) as file_handle:
	for file_line in file_handle:
		row_counter = row_counter + 1
		file_line = file_line.strip("\n")
		if row_counter <= args.header_rows:
			print(file_line)
		else:
			file_line_array = file_line.split("\t")	
			chrom_name = file_line_array[args.chrom_col - 1] # python arrays start at 0
			if chrom_name[0:5] == 'group' :
				chrom_name = 'chr' + chrom_name[5:999]
				file_line_array[args.chrom_col - 1] = chrom_name
			if chrom_name[0:8] == 'scaffold' :
				# if the chrom name did not begin with group, we will use the conversion dataframe
				# new_chrom_name = (convert_df.loc[convert_df['Scaffold']==chrom_name]['Chr'].values)[0] ## use column names passed as args below
				new_chrom_name = (convert_df.loc[convert_df[args.cf_frm_name_col]==chrom_name][args.cf_to_name_col].values)[0]
				# read the offset from the conversion dataframe
				# offset = (convert_df.loc[convert_df['Scaffold']==chrom_name]['Start'].values)[0] ## use column names passed as args below
				offset = (convert_df.loc[convert_df[args.cf_frm_name_col]==chrom_name][args.cf_start_col].values)[0]
				offset = offset - 1 # subtract one from the offset
				# determine the new_start value
				start = file_line_array[args.start_col - 1]
				new_start = int(start) + offset
				# determine the new_end value
				end = file_line_array[args.end_col - 1]
				new_end = int(end) + offset
				# put the new data back in the file_line_array
				file_line_array[args.chrom_col - 1] = new_chrom_name
				file_line_array[args.start_col - 1] = str(new_start)
				file_line_array[args.end_col - 1] = str(new_end)
				
			# print(chrom_name[0:5]) ## this was just used when testing
			out_string = "\t".join(file_line_array)
			print(out_string)
