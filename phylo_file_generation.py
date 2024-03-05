import pandas as pd
from Bio import SeqIO
import shutil
import os
from datetime import datetime
import subprocess
import rpy2.robjects as robjects

# ---------------------------------------------------------------------------------
# INPUT DATA FILTERING CRITERIA
# ---------------------------------------------------------------------------------

# Define criteria for selecting sample data
sample_filters = {
    "Year": [2021],  # Only include samples from these years
    "Sample origin": "Human"  # Only include samples originating from humans
    # Expand with additional criteria as necessary
}

# Select the year to use as a reference for analysis
analysis_reference_year = 2019

# ---------------------------------------------------------------------------------
# DATA PROCESSING FUNCTIONS
# ---------------------------------------------------------------------------------

def filter_dataframe(df, filters):
    """
    Filters a DataFrame based on specified criteria for each column.

    Parameters:
    - df: pandas.DataFrame to filter
    - filters: dict where keys are column names and values are the required entries
               (either a single value or a list of values)

    Returns:
    - A filtered pandas.DataFrame
    """
    
    filtered_df = df.copy()
    for column, criteria in filters.items():
        if criteria is not None:
            if isinstance(criteria, list):
                # If the criteria is a list, filter to include any of the values in the list
                filtered_df = filtered_df[filtered_df[column].isin(criteria)]
            else:
                # If the criteria is a single value, filter based on that value
                filtered_df = filtered_df[filtered_df[column] == criteria]
    return filtered_df


def create_sample_dates_file(csv_file_path, output_folder, output_file_name):
    """
    Creates a file with Sample ID and Tip.Dates from a CSV, tab-separated.

    Parameters:
    - csv_file_path: Path to the input CSV file.
    - output_folder: Folder where the output file will be saved.
    - output_file_name: Name of the output file.
    """
    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Load the DataFrame from CSV
    df = pd.read_csv(csv_file_path)

    # Extract 'Sample ID' and 'Tip.Dates' columns
    sample_dates_df = df[['Sample ID', 'Tip.Dates']]

    # Specify the file path for the output file
    output_file_path = os.path.join(output_folder, output_file_name)

    # Write the data to the file, tab-separated, without index and header
    sample_dates_df.to_csv(output_file_path, sep='\t', index=False, header=False)
    

def create_districts_file(csv_file_path, output_folder, output_file_name):
    """
    Reads a CSV file, extracts Sample ID and District, and saves to a new text file.
    Adds custom column names 'traits' and 'districts'.

    Parameters:
    - csv_file_path: Full path to the input CSV file.
    - output_folder: Folder where the output file will be saved.
    - output_file_name: Name for the output text file.
    """
    # Load the DataFrame from the specified CSV file
    df = pd.read_csv(csv_file_path)

    # Select and rename the required columns
    districts_df = df[['Sample ID', 'District']].rename(columns={'Sample ID': 'traits', 'District': 'districts'})

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Construct the full path for the output file
    output_path = os.path.join(output_folder, output_file_name)

    # Save the DataFrame to a text file, tab-separated, including the header
    districts_df.to_csv(output_path, sep='\t', index=False, header=True)


### ---------------------------------------------------------------------------------
### CREATE PROJECT FOLDER
### ---------------------------------------------------------------------------------
# Generate a folder name with a timestamp
timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
project_folder = f"projects/data_{timestamp}"
os.makedirs(project_folder, exist_ok=True)


### ---------------------------------------------------------------------------------
### FILTERED METADATA w/ README
### ---------------------------------------------------------------------------------

# Load your DataFrame (replace with the correct path to your CSV file)
all_samples_df = pd.read_csv("data/D1-4_NDA_SM_NUS.metadata.csv")

# Save the README file with the filters description
readme_path = os.path.join(project_folder, "DESCRIPTION.md")
with open(readme_path, 'w') as readme_file:
    readme_file.write("# Data Filters Applied\n\n")
    readme_file.write("This folder contains data filtered with the following criteria:\n")
    for key, value in sample_filters.items():
        readme_file.write(f"- **{key}**: {value}\n")

# Filter the DataFrame
filtered_sample_df = filter_dataframe(all_samples_df, sample_filters)


### ---------------------------------------------------------------------------------
### FASTA
### ---------------------------------------------------------------------------------

# Assuming `filtered_sample_df` is your DataFrame and contains the `Sample ID` column
sample_ids = filtered_sample_df['Sample ID'].tolist()

# Path to your combined FASTA file
combined_fasta_path = "data/D1-4_NDA_SM_NUS.fasta"

# Assuming outfile is your FASTA file name
filtered_fasta_path = os.path.join(project_folder, "sample_data.fasta")

# Filter and write the sequences to the new FASTA file
with open(filtered_fasta_path, 'w') as outfile:
    for record in SeqIO.parse(combined_fasta_path, "fasta"):
        # Check if the record.id is in your list of Sample IDs
        if record.id in sample_ids:
            SeqIO.write(record, outfile, "fasta")
            

### ---------------------------------------------------------------------------------
### Align sequences to reference
### ---------------------------------------------------------------------------------

# Find the most common value (mode) in the 'Serotype' column
most_common_serotype = filtered_sample_df['Serotype'].mode()[0]

references_metadata = pd.read_csv('data/references/references.csv')

# Find the row in references_metadata that matches the conditions
matching_row = references_metadata[
    (references_metadata['Year'] == analysis_reference_year) & 
    (references_metadata['Serotype'] == most_common_serotype)
]

print(matching_row)

# Since there might be multiple rows matching the condition, we'll take the first match
# If you expect only one match or want the first match, use .iloc[0:1] to get it as a DataFrame
matching_row = matching_row.iloc[0:1]

# Concatenate the matching row with filtered_sample_df, making sure the matching row is first
# Reset the index to account for the new row order and drop the old index
filtered_sample_df = pd.concat([matching_row, filtered_sample_df], ignore_index=True)

# Specify the file path including the folder name
metadata_path = os.path.join(project_folder, "dengue_samples.metadata.csv")

# Assuming filtered_sample_df is your filtered DataFrame
filtered_sample_df.to_csv(metadata_path, index=False)


### ---------------------------------------------------------------------------------

# Assuming matching_row is defined and contains a 'Sample ID' column
sample_id = matching_row['Sample ID'].iloc[0]  # Get the sample ID value

# Path to your FASTA file
fasta_file = 'data/references/references.fasta'

# Initialize variable to store the matching sequence
matching_sequence = None

# Open and iterate through the FASTA file to find the matching sequence
with open(fasta_file, 'r') as fasta:
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id == sample_id:
            matching_sequence = record
            break  # Exit the loop once the matching sequence is found


reference_fasta_path = os.path.join(project_folder, f'reference_{analysis_reference_year}_{most_common_serotype}.fasta')


# Save the matching sequence to the new FASTA file
with open(reference_fasta_path, 'w') as output_file:
    SeqIO.write(matching_sequence, output_file, 'fasta')


### ---------------------------------------------------------------------------------

# Run the Bash script
result = subprocess.run(['bash', "./MAPLE_Tree.sh", "sample_data.fasta", 
                         f'reference_{analysis_reference_year}_{most_common_serotype}.fasta', 
                         project_folder], capture_output=True, text=True)


### ---------------------------------------------------------------------------------
# MAKE BEAUTI FOLDER
    
# Create the folder
os.makedirs(f"{project_folder}/beauti", exist_ok=True)

### ---------------------------------------------------------------------------------
### sampleDates.txt, sampleDistricts.txt
### ---------------------------------------------------------------------------------

create_sample_dates_file(metadata_path, f"{project_folder}/beauti", "sampleDates.txt")
create_districts_file(metadata_path, f"{project_folder}/beauti", "districts.txt")

# Define the source and destination paths
source_path = f"{project_folder}/sample_data.MAPLE/sample_data.aligned.fasta"
destination_path = f"{project_folder}/beauti/input_sequences.fasta"

# Ensure the destination folder exists, create if it doesn't
destination_folder = os.path.dirname(destination_path)
os.makedirs(destination_folder, exist_ok=True)

# Copy the file and rename it
shutil.copy(source_path, destination_path)


## ---------------------------
## TFPSCANNER
## ---------------------------

## TURN METADATA TO TRAITS.CSV
## MAKE FOLDER FOR DIFFERENT REFERENCES, INCLUDE METADATA
## RUN 
    
# # Define R script and run it
# r_script = """

# # =setwd("/Users/conor/NUS/dengue/phylodynamics_dengue_singapore")

# # metadata <- read.csv("projects/data_20240301-110654/traits.csv")
# # tree <- read.tree("projects/data_20240301-110654/filtered.MAPLE/filtered.MAPLE_output_tree.tree")

# # tfpscanner::tfpscan(tre = tree, amd = metadata, min_descendants = 5)

# # tfpscanner::create_browser_data(
# #   e0 = "tfpscan-2024-03-01/scanner-env-2024-03-01.rds",
# #   output_dir = "tfpbrowser_files"
# # )
# """

# # Run R script
# robjects.r(r_script)