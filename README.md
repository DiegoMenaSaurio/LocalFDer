# LocalFDer

LocalFDRer is a Python script that calculates the local FDR from a subpopulation of PSMs with a 
specific type of post-translational modification using PSMs.txt files from Proteome Discoverer 2.5 as input

**Syntax**:LocalFDRer.py -cConfig_LocalFDRer.txt

The necessary parameters to execute LocalFDRer must be indicated in the "Config_LocalFDRer.txt" configuration file:
* **infile** is the path to the directory containing the PSMs.txt files;
* **outfile** is the path to the output directory containing the PSMs remaining after local FDR filtering;
* **finalname** is the output file's name;
* **aminoacids_list** indicates which residues (separated by comma) must be considered for the calculation of the local FDR;
* **modification_to_count** is the Unimod modification on amino_acids to be considered for the calculation of the local FDR;
* **amino_acids_to_exclude** indicates which residues (separated by comma) should not be considered for the calculation of the local FDR;
* **modification_to_exclude** is the Unimod modification on amino_acids_to_exclude to be excluded from the calculation of the local FDR;
* **FDR** is the local FDR threshold to be considered for populating outfile. Together with outfile, 

LocalFDRer generates a Summary file indicating the number of PSMs obtained using local and global FDR as well as the percentage local PSMs over global PSMs.
