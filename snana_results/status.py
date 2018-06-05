#!/usr/bin/env python

from glob import glob

file_list = sorted(glob('*.log'))
max_name_length = max([len(name) for name in file_list]) - 4

for file_name in file_list:
     with open(file_name, 'r') as ofile:
         end_status = ofile.readlines()[-1].strip()
         
         out_msg = file_name.strip('.log').ljust(max_name_length)
         if end_status == 'ENDING PROGRAM GRACEFULLY.':
             out_msg += ' : Ended Gracefully'
         
         else:
             out_msg += ' : !! Early Exit'
             
         print(out_msg)

print()
