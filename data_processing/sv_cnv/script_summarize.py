'''import os
import csv
import sys

# Get the directory from command line arguments
directory = sys.argv[1]

# Iterate over all files in the specified directory
for filename in os.listdir(directory):
    if filename.endswith('_plink_out'):
        with open(os.path.join(directory, filename), 'r') as f:
            reader = csv.reader(f, delimiter=',')
            new_rows = []

            # Skip header and add new header names
            header = next(reader)
            header.extend(['het_count', 'hom_count'])
            new_rows.append(header)

            for row in reader:
                # Calculate het_count and hom_count, skipping the first column
                het_count = row[1:].count('1')
                hom_count = row[1:].count('2')

                row.extend([str(het_count), str(hom_count)])
                new_rows.append(row)

        # Write the processed rows to the new file
        new_filename = filename.replace('_plink_out', '_plink_out_processed')
        with open(os.path.join(directory, new_filename), 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(new_rows)'''

import os
import csv
import sys

# Get the directory from command line arguments
directory = sys.argv[1]

# Iterate over all files in the specified directory
for filename in os.listdir(directory):
    if filename.endswith('plink_out'):
        with open(os.path.join(directory, filename), 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader)

            new_rows = []

            # New header names
            header = ['plate_key', 'het_count', 'hom_count']
            new_rows.append(header)

            for row in reader:
                # Calculate het_count and hom_count, skipping the first column
                het_count = row[1:].count('1')
                hom_count = row[1:].count('2')

                # Write only the sample name, het_count, and hom_count to the new row
                new_row = [row[0], str(het_count), str(hom_count)]
                new_rows.append(new_row)

        # Write the processed rows to the new file
        new_filename = filename.replace('plink_out', 'plink_out_processed')
        with open(os.path.join(directory, new_filename), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(new_rows)

