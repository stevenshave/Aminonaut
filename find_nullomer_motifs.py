
"""
Identify motifs which account for nullomer peptides

Reading in files containing counts of unique peptides at certain lengths,
generate motifs (20 AAs + any AA) and count how many times the motif
accounts for a nullomer peptide.  This allows identification of sequences that
life does not like to make.
"""

__author__ = "Steven Shave"
__version__ = "1.0.0"
__license__ = "MIT"

import argparse, re, gzip
from itertools import product


def fine_nullomer_motifs(input_filename:str, output_filename, pattern_length):
    """ 
    Function to explore all possible motifs of lenght pattern_length and count
    how well they define nullomers present in input_filename.  Writes csv to
    output_filename.
    """
    
    # Valid amino acids plus . to represent any AA.  Dot is regex for any char.
    aas="ARNDCEQGHILKMFPSTWYV."
   
    # Enumerate all combinations of chars in aas. Patterns is used to make sure
    #  we have not seen it before
    patterns=set()
    unique_patterns=[x for x in product(*[aas]*pattern_length) if tuple(x[::-1]) not in patterns and not patterns.add(tuple(x))]
    print("Unique patterns", unique_patterns)
    
    # Compile to regular expressions all combinations of chars 
    regex_queries=[re.compile("".join(m)) for m in unique_patterns]

    # Occurences dictionary holds number of pattern matches
    occurences={}
    print("len queries=", len(regex_queries))
    nullomer_counter=0
    
    input_file=None
    if input_filename[-3:]==".gz":
        input_file=gzip.open(input_filename, 'rt')
    else:
        input_file=open(input_filename)

    for line_it, line in enumerate(input_file.readlines()):
        if line_it%1000==0: # Progress counter
            print(line_it)
        if line.find(",")==-1:continue # If line does not contain a comma, skip
        if int(line.split(",")[1])>0:continue # If not a nullomer, skip
        nullomer_counter+=1
        for regex_query in regex_queries:
            regex_result=regex_query.findall(line)+regex_query.findall(line[::-1])
            if len(regex_result)>0: # If regex matches, add to dict or increment
                if regex_query.pattern in occurences.keys():
                    occurences[regex_query.pattern]+=len(regex_result)
                else:
                    occurences[regex_query.pattern]=len(regex_result)
    input_file.close()
    output_file=None
    # Write the output file, containing regex, number of matches.
    if output_filename[-3:]==".gz": # Writing compressed gz file
        output_file=gzip.open(output_filename,"wb")
        output_file.write(f"{'Motif,':>10}{'Count,':>10}{'%Match,':>10}\n".encode())
        for l in sorted(occurences.items(), key=lambda x: x[1], reverse=True):
            output_file.write(f"{l[0]+',':>10}{str(l[1])+',':>10}{(100*l[1]/nullomer_counter):>8.3f}%\n".encode())
        output_file.close()
    else:
        output_file=open(output_filename, "w")
        output_file.write(f"{'Motif,':>10}{'Count,':>10}{'%Match,':>10}\n")
        for l in sorted(occurences.items(), key=lambda x: x[1], reverse=True):
            output_file.write(f"{l[0]+',':>10}{str(l[1])+',':>10}{(100*l[1]/nullomer_counter):>8.3f}%\n")
        output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("input_filename", help="File containing peptides and counts comma separated")
    parser.add_argument("output_filename",
                        help="File to output motif hit counts to")
    parser.add_argument("motif_length",
                        help="File containing sequences", type=int)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    fine_nullomer_motifs(args.input_filename, args.output_filename,
                       args.motif_length)