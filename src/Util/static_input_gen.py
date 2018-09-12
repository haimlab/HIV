""" 
generate csvs separated by clade and region for static data
primarily used in profile centroid p-value calculations
"""

import csv
from os.path import join

in_file_name = "constants/Copy of 8.29.18 2G12-'07'-15 - Clades_Positions_Geo.csv"
out_dir = ''
amino_acids = ['Z', 'N', 'T', 'S', 'D', 'E', 'K', 'R', 'H', 'Y', 'Q', 'I', 'L', 'V', 'A', 'C', 'F', 'G', 'M', 'P', 'W']

with open(in_file_name) as file:
    reader = csv.reader(file)
    buffer = []
    for row in reader:
        buffer.append(row)
        if len(buffer) == 25:
            rc = buffer[0][0]
            positions = [buffer[1][i] for i in range(0, 13, 3)]
            fns = [rc + '_' + p + '.csv' for p in positions]
            
            i = 1
            for fn in fns:
                with open(join(out_dir, fn), 'w') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    for aa, cur_row in zip(amino_acids, buffer[2:23]):
                        writer.writerow([aa, cur_row[i]])
                i += 3
            
            buffer = []
                        
        
        
        
        
        
        
        
        #if i % 25 == 0:
            #clade_region = row[0]
        #elif i % 25 == 1:
            #positions = [row[m] for m in range(0, 13, 3)]
        #elif i % 25 != 24:
            #distr = []
            #for k in range(0, 21):
                #amino_acid = row[0]
                #distr += [row[m] for m in range(1, 14, 3)]
        #else:
            #file_names = [cr + p + '.csv' for cr, p in zip(clade_region, positions)]
            #f_ind = 0
            #for fn, distr in zip(file_names, distr):
                #with open(fn, 'w') as f:
                    
                    
                    
        