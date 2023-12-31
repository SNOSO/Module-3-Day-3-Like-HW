# Module-3-Day-3-Like-HW

README

Given a full network(STRING1.txt) and a gene-locus (Input.GMT.txt) file, compute the score of each gene at the locus and visualize the networks from the generated data frames. 

INSTALLATION AND DEPENDENCIES: This program was written using Python 3.11.2 and requires the following packages: 
pandas
numpy
random

installed using pip. For example: pip install <package>

USAGE: 

Open the file in a code text editor of your choice and input the path to your files in the designated file_name inputs where:
string = your full network
FA_genes = genes of your target disease

OPTIONS:

-input_full_network name of the full netowrk file where the format follows three columns, column 1 is Gene1, column 2 is Gene2 and column 3 is the weight of the edge between these two genes[STRING.txt]

-input_gene_locus name of the gene-locus GMT file where the first column is the locus name, the rest of the rows are the genes associated with the locus[Input.gmt.txt]

-num_net specify number of random prix fixe network solutions to generate [5000]

-out_file1 specify the name of the FA genes dataframe to use as the network file in Cytoscape[specifed_output_network_name.txt]

-out_file2 specify the name of the genes, scores and locus to use as a node table in Cytoscape [specifed_output_name.txt]

RUN:
To run this program in command line type:

$ python3 Mod3Day3LikeHW.py 

in the same folder as the python file and the input files. Unless specified, the output file will be written to the same directory.
