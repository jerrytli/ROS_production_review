## Part 1: Calculate global full-length NADPH oxidase homologs

 1. Download conda and blast.
    
    We can use conda to help download software needed to find the orthologs we desire and processes them for phylogenetic analsyes.
 
    ```bash
    # go to this website https://www.anaconda.com/distribution/
    # select the 3rd distribution if you don't already use anaconda 
    bash Anaconda-latest-Linux-x86_64.sh


    # run so your computer can find conda with your path
    #for anaconda 2 :

    export PATH=~/anaconda2/bin:$PATH

    #for anaconda 3 :

    export PATH=~/anaconda3/bin:$PATH


    if you need to initalize conda
    #write the command $conda init {shell script - aka bash}
    #then restart ternminal, it should look like this: 
    #(base) danimstevens@danimstevens-MS-7B93:

    #set up other conda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```
   
    After conda is initalized, we can use it to stabaily download a variety of bioinformatic software. 
    Search this website to see what they have available: https://anaconda.org/bioconda
   
    ```bash
     conda install -c bioconda blast 
     ```

 2. Make blast database of genes you want to search for:
 
    First create a text (fasta) file of the genes we want to blast.
    ```bash
    # to make blast db
    makeblastdb -in NADPH_oxidase_homologs.fasta -parse_seqids -dbtype 'prot' -out NADPH-oxidase-homologs-db
    ```
  3. Pull out the c-terminus domain of NADPH oxidase homologs using the C-terminus AtRBOHD as a query.
        ```bash
        blastp -query NADPH-oxidase-cTerm-reference.fasta -db ./BLAST_database/NADPH-oxidase-homologs-db -evalue 1e-8 -outfmt "6 qseqid sseqid pident evalue len qstart qend sseq" -out blast-hits-c-terminus.txt
        ```
        
       I increased the E-value cutoff to limit the number of incomplete or low quality hits as these could unnesecarily skew the sliding window values.

   4. Processes blast outputs into fasta files for alignment: </br>
   
       ```
       I always manually assess the blast results and see if they are reasonable (if not, rerun and alter options as 
       necessary and pertinent). To processes these hits, the .txt file will be converted to a fasta file. Open the 
       convert_blast2fasta_script.py python script to processes the text files into fasta files. Currently this script 
       is not command line friendly and has to run line-by-line. However, it is stable in rstudio/rserver.
       ```
    
   5. Process the fasta file blast hits to determine similarity values.
    
      ```
      I courrently working on editing the code so everything will run on R consel but for now, just open the script 
      and run line-by-line. It should be noted that collecting simiarity values in a scanning-window appraoch takes 
      a while to run (1-2 hours). I am working on reshashing the logic so it runs more efficently in the future, but 
      for now, I added a progress bar for easier progress monitoring.
      ```
   
## Part 2: Calculate similarity across residues of C-terminus AtRBOHD
   
   5. Process the fasta file blast hits to determine similarity values. 
   
     (Same script as in #5 from Part 1).

