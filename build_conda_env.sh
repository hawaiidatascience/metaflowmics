conda create -n cmaiki -y python==3.7.3
conda activate cmaiki
conda install -y nextflow==19.04.1
conda install -y biopython==1.73 pandas==0.24.2
conda install -y r-base==3.5.1
conda install -y r-ggplot2 bioconductor-dada2 r-stringr r-seqinr r-devtools r-dplyr
Rscript -e 'devtools::install_github("tobiasgf/lulu")'
conda install -y mothur==1.41.3 vsearch==2.13.3
