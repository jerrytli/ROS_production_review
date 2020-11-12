## Part 1: Calculate global similarity of full-length NADPH oxidase homologs

1. Download anaconda and blast.
    
    We can use conda to help download software needed to find the orthologs we desire and processes them for phylogenetic analsyes.
 
    ```bash
    # go to this website https://www.anaconda.com/distribution/
    # select the 3rd distribution if you don't already use anaconda 
    # my computer is a linux, so I have downloaded the linux distribution and can use bash to intall it
    bash Anaconda-latest-Linux-x86_64.sh


    # run so your computer can find conda with your path
    # for anaconda 2 :

    export PATH=~/anaconda2/bin:$PATH

    # for anaconda 3 :

    export PATH=~/anaconda3/bin:$PATH


    # if you need to initalize conda
    # write the command $conda init {shell script - aka bash}
    # then restart ternminal, it should look like this: 
    (base) danimstevens@danimstevens-MS-7B93:

    # set up other conda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```
   
    After conda is initalized, we can use it to download a variety of bioinformatic software. 
    Search this website to see what they have available: https://anaconda.org/bioconda
    Note: Sometimes the package versions are not always the most updated, many have developer versions on github
   
    To install the blast library:
    
    ```bash
     conda install -c bioconda blast 
     ```

2. Make blast database of genes you want to search for:
 
    First create a text (fasta) file of the genes we want to blast.
    
    ```bash
    # to make blast db
    makeblastdb -in NADPH_oxidase_homologs.fasta -parse_seqids -dbtype 'prot' -out NADPH-oxidase-homologs-db
    ```
    
3. Pull out the c-terminus and n-terminus domain of NADPH oxidase homologs using the N-terminus and C-terminus AtRBOHD as a query.
  
    ```bash
    # N-terminus hits
    blastp -query NADPH-oxidase-nTerm-reference.fasta -db ./BLAST_database/NADPH-oxidase-homologs-db -evalue 1e-10 -outfmt "6 qseqid sseqid pident evalue len qstart qend sseq" -out blast-hits-n-terminus.txt
      
    # C-terminus hits
    blastp -query NADPH-oxidase-cTerm-reference.fasta -db ./BLAST_database/NADPH-oxidase-homologs-db -evalue 1e-10 -outfmt "6 qseqid sseqid pident evalue len qstart qend sseq" -out blast-hits-c-terminus.txt
    ```
        
    I increased the E-value cutoff to limit the number of low quality hits as these could unnesecarily skew the sliding window values.

4. Processes blast outputs into fasta files for alignment: </br>
   
    ```
    I always manually assess the blast results and see if they are reasonable (if not, rerun and alter options as 
    necessary). To processes these hits, the .txt file will be converted to a fasta file. Open the 
    convert_blast2fasta_script.py python script to processes the text files into fasta files. Currently this script 
    is not command line friendly and has to run line-by-line. However, it is stable in rstudio/rserver.
    ```
    
5. Process the fasta file blast hits to determine similarity values.
    
    ```
    This R script (Process_NADPH_homologs.R) is also ran line-by-line. But there is some additionally documentation (comments) 
    that should make it relatively easy to follow.
    It should be noted that collecting simiarity values in a scanning-window appraoch takes 
    a while to run (1-2 hours). I added a progress bar for easier progress monitoring.
    ```
   
## Part 2: Calculate similarity across residues of the N- and C-terminus AtRBOHD
   
1. Process the fasta file blast hits to determine similarity values. 
   
     (Same script as in #5 from Part 1).

