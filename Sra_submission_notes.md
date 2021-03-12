Notes on Uploading Data to NCBI SRA
-------------------------------------------------
NCBI offers an online wizard to submit sequences for publication, this document offers some notes to clarify the process.

# FTP preload for large files
Navigate to the [Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/sra/) and open the FTP upload tab

1. On the computer the sequences are, Open a terminal in the sequence directory. 

2. Run the command ftp ftp-private.ncbi.nlm.nih.gov 
    * Username: subftp
    * Password : see NCBI website

3. Change directory (cd) into the dir listed on the NCBI website
    * cd uploads/<$STRING_SPECIFIC_TO_ACCT>

4. Create a new directory to store sequences from this project
    * mkdir $PROJECT

5. Change directory (cd) into that new folder
    * cd $PROJECT

6. Copy sequences into that new directory
    * put $FILES

# Ordering of steps
* FTP preload
* Create BioProject

* Create BioSample 

* SRA Submission Portal

