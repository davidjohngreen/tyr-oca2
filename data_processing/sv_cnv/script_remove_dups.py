import os
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python3 script.py <file1> <file2>")
    sys.exit(1)

# Assign input arguments to variables to get input and output filenames
file1 = sys.argv[1]
file2 = sys.argv[2]


# Keep track of which lines we've seen before
seen_lines = set()

# Open the input and output files
with open(file1, 'r') as f_in, open(file2, 'w+') as f_out:
    # Loop through the lines in the input file
    for line in f_in:
        if line.startswith('#'):
            # If the line starts with a "#" character, write it to the output file
            f_out.write(line)
        else:
            # Otherwise, split the line into fields
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            id = fields[2]

            # Check if we've already seen this line
            line_key = "{}:{}:{}".format(chrom, pos, id)
            if line_key not in seen_lines:
                # If we haven't seen this line before, write it to the output file
                f_out.write(line)
                seen_lines.add(line_key)
