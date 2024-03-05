# Phylogeographic Insights into the Dengue Outbreak in Singapore: Tracing Origins and Predicting Spread through High-Throughput Sequencing

This script generates files used for phylogenetic/dynamic/geographic analysis using the BEAST software. The data is omitted, but is available on my OneDrive (can place elsewhere if there is a particular preference). 

## Data

The sequencing data is sourced from NEA Singapore. This includes samples for the four dengue serotypes, between the years 2013 and 2022, with its corresponding metadata that includes location, receipt-date, genotype, etc. In the `data` folder (available on my OneDrive), there are the following folders:

- A folder for each serotype, containing the fasta file of all sequences of said serotype, and a csv metadata file
- Case distribution information - cases per year, and per month
- Reference sequences (with corresponding metadata), one per year (2012-2022) and per serotype

as well as the sequence fasta file for all samples and metadata.

Place this data folder in the cloned repository for use.

## Other files/folders

- The MAPLE software used in the relatively fast creation of a tree given a set of sequences and a reference
- A project folder, where all analyses end up in their respective folders (projects will have a time-stamped name)
- `phylo_file_generation.py`, a Python script used to generate the following files/folders, based on the specified filters applied to the entire dataset - e.g. range of years, district origin, etc.:
    - A description md file which lists the filters applied
    - Filtered sequence fasta file (with metadata)
    - Reference file selected (choose the year from 2012 to 2022)
    - Folder containing generated MAPLE files (aligned sequences to reference, tree-file, substitution matrix, and input MAPLE file for tree generation)
    - BeaUTi input files - district traits, time of samples, and the input fasta file

## Running the programme

At the top of `phylo_file_generation.py`, you need to define the filter criteria. This can be the required years, sample origin, etc. Have a look at the headers of the metadata, and select the required values. If multiple required, include as a list.

Also select the year you want the reference from. Generally, we take the reference to be the year before the earliest selected sample. The programme will automatically choose the serotype corresponding to the most prevalent serotype in the filtered dataset. 


## Future extensions

If it is possible to run BeaUTi without having to use the GUI, then the natural extension would be to have it run in the script, after specifying the prior, substitution model, MCMC parameters, etc.

I also want to include the analysis performed by `tfpscanner`. This requires the tree file, and a traits file (Sample ID, district location, receipt-date). I don't have this working yet - hopefully soon! 
