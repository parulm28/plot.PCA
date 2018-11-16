
*** The R Code was developed by Parul Mittal.

The package can be used to carry out Principal Component Analysis and returns a PCA 
plot. The package was developed to identify if the time points at which the samples 
are collected for an analysis differentiates on a PCA plot. The package identifies 
how a given set of features varies at different time points. For PCA, the package 
utilizes prcomp package, so the PCA calculation will be done by a singular value 
decomposition of the centered and scaled data matrix.

Input:

File 1: on which PCA analysis is to be carried out. The file should be csv file with 
the first row as Sample names. The second column of the file should contain unique IDs 
with a header "symbol". The third to the last column of this file should contain (numeric) values corresponding to that feature. For example:

,symbol,sample1,sample2,sample3
1,gene1,12.2222,188.876,23.8777
2,gene2,32.4226,128.172,73.7657

The values 12.2222, 188.876, 23.8777 are relative gene expression for gene1.
The values 32.4226, 128.172, 73.7657 are relative gene expression for gene2.

It can have any number of rows(genes) and columns(samples).

File 2: is the Metadata file. It should also be a csv file with the first row as headers.
The first column should contain the Sample names as given as the headers in File 1. 
The header of this column should be "sIdx". One column of this should contain Time column. 
This column is used to distinguish the data.

Example of Metadata file:

sIdx,Time
sample1,2
sample2,4
sample3,8

If there is any junk/missing value for any gene, that gene will be eliminated from the 
analysis. If more than 80% genes are removed because of this, the program will be 
terminated.

The package requires ggplot2 for plotting. It is required that you install the package
in R by using the command:

install.packages("ggplot2")

The package also uses prcomp function which is available in the stats package.

How to use:

Rscript script.R <file1.csv> <file2.csv>



