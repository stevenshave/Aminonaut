# NullomerPeptides
Software to identify and analyse nullomer peptides from UniProt (Swiss-Prot) XML data

This software is an output of the Leverhulme funded project:
_"Understanding the biological ramifications of a ‘forbidden’ peptide nullomer."_

## The programs
To generate counts of observed peptides within nature, identify the nullomer peptides and then perform motif searching, first obtain the most recent UniProt Swiss-Prot release as XML. Then in succession and with the usage bellow, run the count_peptides and find_nullomer_motifs program. The software tools contained in this package support compressed input and output, inferred from the filename.  Files ending in .gz are considered gzip compressed.

### Count_peptides.py
This program reads in an XML file such as UniProt (Swiss-Prot), extracting protein sequences between \<sequence\> tags.  Once found, and operating on a certain length peptide, the number of unique peptides of that length are obtained from the file.

count_peptides.py was used to identify peptide nullomers in the _"Understanding the biological ramifications of a ‘forbidden’ peptide nullomer"_ project.

#### Usage
```
python count_peptides.py

usage: count_peptides.py [-h] [-c OUTPUT_CUTOFF] [--version] uniprot_input_file output_filename peptide_lengths
count_peptides.py: error: the following arguments are required: 
uniprot_input_file, output_filename, peptide_lengths
```
#### Extended usage
```
python .\count_peptides.py -h
usage: count_peptides.py [-h] [-c OUTPUT_CUTOFF] [--version] uniprot_input_file output_filename peptide_lengths

positional arguments:
  uniprot_input_file    File containing sequences
  output_filename       File to output sequence counts to
  peptide_lengths       Length of peptides

optional arguments:
  -h, --help            show this help message and exit
  -c OUTPUT_CUTOFF, --output_cutoff OUTPUT_CUTOFF
                        Dont output sequences apearing more than this many times
  --version             show program's version number and exit
```

To generate 2,3,4,5, and 6-mer peptide counts, the following commands may be used:
> python count_peptides.py uniprot_swissprot.xml.gz petides-2mers.txt 2

> python count_peptides.py uniprot_swissprot.xml.gz petides-3mers.txt 3

> python count_peptides.py uniprot_swissprot.xml.gz petides-4mers.txt 4

> python count_peptides.py uniprot_swissprot.xml.gz petides-5mers.txt.gz 5

> python count_peptides.py uniprot_swissprot.xml.gz petides-6mers.txt.gz 6

Note we specify a file extension ending in .gz for the 5 and 6-mers. This is because of the large number of unique peptides possible in 5 and 6-mer peptides (20^5, and 20^6).



### find_nullomer_motifs.py and find_nullomer_motifs_forward_and_backwards.py
These program read in output from the previous count_peptides.py program, essentailly CSV files with the first column containing the unique peptide sequence, and followed by another column with the number of times the peptide was found.  With this read in, it constructs set length motifs and queries how many nullomers this motif covers.  For example, with the identification of CQWW, we may theorise that the motiff C..W where dot is any amino acid is enriched within the 5-mers.  We may use find_nullomer_motifs to enumerate all possible motifs and query the previously generated datasets, counting how many nullomers match the queried motif. The program find_nullomer_motifs matches motifs in only the forwards direction, so C..W would match CQWW, and not WWQC.  The program find_nullomer_motifs_forward_and_backwards matches in both directions, so the morif C..W would match CQWW and WWQC.

#### Usage
```
python find_nullomer_motifs.py
usage: find_nullomer_motifs.py [-h] [--version] input_filename output_filename motif_length
find_nullomer_motifs.py: error: the following arguments are required: input_filename, output_filename, motif_length
```
#### Extended usage
```
python find_nullomer_motifs.py -h
usage: find_nullomer_motifs.py [-h] [--version] input_filename output_filename motif_length

positional arguments:
  input_filename   File containing peptides and counts comma separated
  output_filename  File to output motif hit counts to
  motif_length     File containing sequences

optional arguments:
  -h, --help       show this help message and exit
  --version        show program's version number and exit
```

To generate generate motifs and their nullomer coverage, we cann use the following:
To identify nullomer motifs of length 4 within the 5-mer peptides:
> python find_nullomer_motifs petides-5mers.txt.gz motifs-4mersIn5mers.txt 4

To identify nullomer motifs of length 3 within the 6-mer peptides:
> python find_nullomer_motifs petides-6mers.txt.gz motifs-3mersIn6mers.txt 3


## Requirements
* Python (>=3.6)
* Numpy (>=1.15.0)





