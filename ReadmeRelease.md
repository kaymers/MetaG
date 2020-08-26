# Configuring

1. Download MetaG_with_backend.rar and extract it.

2. Install python dependancies of the repository

````
(in MetaG directory)
pip install -r requirements.txt
````

If you do not have blast configured in your machine follow the below steps.

3. Install BLAST+ executables [manual](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

4. Add blast+ executables insallation directory to *path*.

5. Create a blast database in anywhere in your file system (or configure a standard database like nt, check the NCBI manuals for this).

6. Extract the file rankedlineage.dmp from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip to anywhere in your file system. 

## Configuring the backend

7. Open *./MetaG/MetaG_backend/config.ini* in a text editor

8. Set *blast_db_directory* and *blast_db_name* according to the blast setup in your machine. Look at the sample values already in the file.

9. Set your rankedlineage.dmp path for the field *rankedlineage_path* 

# Running

````
(in MetaG directory)
python app.py
````
