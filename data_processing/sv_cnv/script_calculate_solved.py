import os

# Specify the path where the files are stored
folder_path = "outputs"

# List of genes for which sex-specific processing should be performed
X_linked_genes = ['GRP143', 'FRMD7']

# Create a list to hold the samples to exclude
samples_to_exclude = []

# Loop over all files in the specified directory
for filename in os.listdir(folder_path):
    # Check if this is a .csv file
    if filename.endswith("_plink_out_processed"):
        # Split the filename to get the prefix
        prefix = filename.split("_plink_out_processed")[0]
        print('processing gene... ' + prefix)
        
        # Open the file
        with open(os.path.join(folder_path, filename), 'r') as f:
            lines = f.readlines()  # read the file line by line into a list

        # Skip the first line and split the rest
        data = [line.strip().split() for line in lines[1:]]

        ### PROCESS X-LINKED GENES
        if prefix in X_linked_genes:
            # For each line in data, get the last three elements and assign them to variables
            for line in data:
                sample = line[0]
                het_count = int(line[1])
                hom_count = int(line[2])
                sex = line[3]
                
                # Check if sex is 'male' and het_count > 0
                if sex.lower() == 'male' and het_count > 0:
                    # Check if the sample is already in the list
                    if sample not in samples_to_exclude:
                        # Add the sample to the exclude list
                        samples_to_exclude.append(sample)
                # Check if sex is 'female' and het_count > 1 or hom_count > 0
                if sex.lower() == 'female' and (het_count > 1 or hom_count > 0):
                    # Check if the sample is already in the list
                    if sample not in samples_to_exclude:
                        # Add the sample to the exclude list
                        samples_to_exclude.append(sample)


        elif prefix == 'TYR':
            for line in data:
                # Get elements and assign them to variables
                sample = line[0]
                variant_402 = int(line[1])
                variant_192 = int(line[2])
                promotor = int(line[3])
                het_count = int(line[4])
                hom_count = int(line[5])
                

                # Check the regular conditions for solved via rare variants
                if het_count > 1 or hom_count > 0:
                    if sample not in samples_to_exclude:
                        samples_to_exclude.append(sample)

                # Homozygous haplotypes; remember that for promotor, reference is the effect allele
                if promotor == 0 and variant_402 == 2:
                    if sample not in samples_to_exclude:
                        samples_to_exclude.append(sample)
                
                # Haplotype + a rare variant
                if (promotor == 0 and variant_402 == 1 and het_count > 0) or \
                   (promotor == 1 and variant_402 > 0 and het_count > 0):
                    # Check if the sample is already in the list
                    if sample not in samples_to_exclude:
                        # Add the sample to the exclude list
                        samples_to_exclude.append(sample)


        elif prefix == 'OCA2':
            for line in data:
                # If prefix is not 'TYR' but in main_genes, get the last two elements from fifth and fourth last columns
                sample = line[0]
                het_count = int(line[1])
                hom_count = int(line[2])
                

                try:
                    OCA2_443_count = int(line[3])
                except ValueError:
                    OCA2_443_count = None

                
                if OCA2_443_count == 2:
                    if sample not in samples_to_exclude:
                        samples_to_exclude.append(sample)

                if OCA2_443_count == 1 and het_count > 0:
                    if sample not in samples_to_exclude:
                        samples_to_exclude.append(sample)
                
                # Check if het_count > 1 or hom_count > 0
                if het_count > 1 or hom_count > 0:
                    # Check if the sample is already in the list
                    if sample not in samples_to_exclude:
                        # Add the sample to the exclude list
                        samples_to_exclude.append(sample)

        # This will apply the normal rules to the other autosomal genes besides TYR and OCA2
        else:
            for line in data:
                sample = line[0]
                het_count = int(line[1])
                hom_count = int(line[2])
                
                
                # Check if het_count > 1 or hom_count > 0
                if het_count > 1 or hom_count > 0:
                    # Check if the sample is already in the list
                    if sample not in samples_to_exclude:
                        # Add the sample to the exclude list
                        samples_to_exclude.append(sample)


print('removed due to solved status' + str(len(samples_to_exclude)))

solved_by_rare_and_common = open('solved_rare_common.txt', 'w+')
for sample in samples_to_exclude:
    solved_by_rare_and_common.write(sample + '\t' + 'solved' + '\n')


# SECOND PASS TO CATCH HETS IN TYR AND OCA2
for filename in os.listdir(folder_path):
    if filename.endswith("_plink_out_processed"):
        # Split the filename to get the prefix
        prefix = filename.split("_plink_out_processed")[0]
        if prefix == 'OCA2':
            # Open the file
            with open(os.path.join(folder_path, filename), 'r') as f:
                lines = f.readlines()  # read the file line by line into a list

            # Skip the first line and split the rest
            data = [line.strip().split() for line in lines[1:]]

            for line in data:
                sample = line[0]
                het_count = int(line[1])
                hom_count = int(line[2])
                
                if sample not in samples_to_exclude and het_count == 1:
                    samples_to_exclude.append(sample)


        elif prefix == 'TYR':
            # Open the file
            with open(os.path.join(folder_path, filename), 'r') as f:
                lines = f.readlines()  # read the file line by line into a list

            # Skip the first line and split the rest
            data = [line.strip().split() for line in lines[1:]]

            for line in data:
                sample = line[0]
                het_count = int(line[4])
                hom_count = int(line[5])

                if sample not in samples_to_exclude and het_count == 1:
                    samples_to_exclude.append(sample)



# Print out the samples to exclude

output = open('samples_to_exclude.csv', 'w+')
for sample in samples_to_exclude:
    output.write(sample + '\n')


print('removed after het...' + str(len(samples_to_exclude)))

