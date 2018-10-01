import csv
from os.path import join

in_file_name = "Copy of 9.14.18 2F5-'07'-15 - Clades_Positions_Geo.csv"
out_dir = 'out/'
amino_acids = ['Z', 'N', 'T', 'S', 'D', 'E', 'K', 'R', 'H', 'Y', 'Q', 'I', 'L', 'V', 'A', 'C', 'F', 'G', 'M', 'P', 'W']

with open(in_file_name) as file:
    reader = csv.reader(file)
    buffer = []
    for row in reader:
        buffer.append(row)
        if len(buffer) == 25:
            rc = buffer[0][0]
            if rc[0] == 'B':  # clade B has 6 periods while others have 5
                num_iso = buffer[1][1:7]
                positions = [buffer[2][i] for i in range(0, 29, 8)]
                year_range = buffer[2][1:7]
                fns = [rc.upper() + '_' + p + '.csv' for p in positions]
                
                i = 1
                for fn in fns:
                    with open(join(out_dir, fn), 'w') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow([out_dir[out_dir.rfind('_') + 1:out_dir.rfind('_') + 4]] + year_range)
                        writer.writerow(['num_isolates'] + num_iso)
                        for aa, cur_row in zip(amino_acids, buffer[3:]):
                            writer.writerow([aa, cur_row[i:i+6]])
                    i += 8
                buffer = []
            else:
                num_iso = buffer[1][1:6]
                positions = [buffer[2][i] for i in range(0, 25, 7)]
                year_range = buffer[2][1:6]
                fns = [rc.upper() + '_' + p + '.csv' for p in positions]
                
                i = 1
                for fn in fns:
                    with open(join(out_dir, fn), 'w') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow([out_dir[out_dir.rfind('_') + 1:out_dir.rfind('_') + 4]] + year_range)
                        writer.writerow(['num_isolates'] + num_iso)
                        for aa, cur_row in zip(amino_acids, buffer[3:]):
                            writer.writerow([aa, cur_row[i:i+5]])
                    i += 7
                buffer = []                