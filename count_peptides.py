
"""
Count peptide occurrences from uniprot data

Using a downloaded representation of UniProt (Swiss-Prot) as XML, parse out
all peptides, summing occurence counts of all unique peptides
"""

__author__ = "Steven Shave"
__version__ = "1.0.0"
__license__ = "MIT"

from sys import getsizeof
import argparse
import numpy as np
import gzip
class CodonCounter():
    aa_to_codon_count={
        'A':4,
        'R':6,
        'N':2,
        'D':2,
        'C':2,
        'Q':2,
        'E':2,
        'G':4,
        'H':2,
        'I':3,
        'L':6,
        'K':2,
        'M':1,
        'F':2,
        'P':4,
        'S':6,
        'T':4,
        'W':1,
        'Y':2,
        'V':4,
        '.':61,
    }
    def queryCodonCount(self, motif):
        return [self.aa_to_codon_count[c] for c in motif]
    
def count_peptides(input_file_name, output_file_name, nullomer_length, maximum_count_cutoff:False):
    """
    Count all AA stretches present in input_file_name, write out counts to csv.
    
    Read in a UniProt (Swiss-Prot) XML file, extracting <sequence> records, and counting
    the number of times each possible unique peptide is seen for a given length.
    """

    # Amino acid lookup list, and AA->int and reverse dictionaries for array
    # addressing
    valid_amino_acid_symbols_set = set('FLIMVSPTAYHQENKDCWRG')
    valid_amino_acid_symbols_list = list(valid_amino_acid_symbols_set)
    aa_to_int = {v: valid_amino_acid_symbols_list.index(
        v) for v in valid_amino_acid_symbols_list}
    int_to_aa = {valid_amino_acid_symbols_list.index(
        v): v for v in valid_amino_acid_symbols_list}
    
    # Codon counter - outputs how many unique codons 
    codon_counter=CodonCounter()
    # count the number of lines in the input file
    lines_in_file=None
    if input_file_name[-3:]==".gz":
        lines_in_file = sum(1 for line in gzip.open(input_file_name, "rt"))
    else:
        lines_in_file = sum(1 for line in open(input_file_name))
    print(f"{lines_in_file} lines in {input_file_name}")

    # This is our big n-dimensional array, with each axis representing 
    # a peptide AA position, and of length 20 representing each amino acid
    # A peptide length of 2 would have shape (20, 20), and length 3 would
    # be (20,20,20) etc.
    counts = np.zeros(
        tuple(*[[len(valid_amino_acid_symbols_list)]*nullomer_length]), dtype=int)
    print(
        f"Made lookup array, size is {len(counts)} elements, shape is: {counts.shape}, {getsizeof(counts)/1024/1024:>10.3f} MB")
    
    # Open the file - handling gz compressed if needed
    input_file=None
    if input_file_name[-3:]==".gz":
        input_file=gzip.open(input_file_name, "rt")
    else:
        input_file=open(input_file_name)

    # Read the file
    for line_index, line in enumerate(input_file):
        if line_index % 10000 == 0:
            print(
                f"Progress: {(line_index+1)/lines_in_file*100:>10.2f}%   {line_index:}/{lines_in_file}")
        # Identify sequences
        pos1 = line.find("<sequence")
        if pos1 < 0:
            continue
        pos2 = line.find(">", pos1)
        if pos2 < 0:
            continue
        pos2 += 1
        pos3 = line.find("</sequence", pos2)
        if pos3 < 0:
            continue
        sequence = line[pos2: pos3]
        if len(sequence) < nullomer_length:
            continue

        # Dont include sequences containing invalid characters, such as X,Z,B,J,U, and O which regularly appear.
        illegal_characters = set(sequence)-valid_amino_acid_symbols_set
        if(len(illegal_characters) > 0):
            continue

        # Increment the count for each peptide we see
        for peptide in [sequence[i:i+nullomer_length] for i in range(0, (len(sequence)-nullomer_length)+1)]:
            counts[tuple([aa_to_int[c] for c in peptide])] += 1

    print("Completed readin, now outputting")
    output_file = None
    compressing=False
    header_line="Peptide,OccurrenceRate,PeptideOccurrenceCount\n"
    if output_file_name[-3:]==".gz":
        output_file=gzip.open(output_file_name, "wb")
        compressing=True
        output_file.write(header_line.encode())
    else:
        output_file=open(output_file_name, "w")
        output_file.write(header_line)

    highest_count = np.max(counts)
    overall_count=np.sum(counts)
    print("Overall sum", overall_count)
    print(f"MAX = {highest_count}")
    for current_count in np.unique(counts)[::-1]:
        if current_count<maximum_count_cutoff: break
        print("Outputting, ", current_count)
        indexes_of_current_count = np.argwhere(counts == current_count)
        for counts_index in indexes_of_current_count:
            peptide = "".join([int_to_aa[i] for i in counts_index])
            codon_counts=codon_counter.queryCodonCount(peptide)
            codon_chance=np.prod([x/61 for x in codon_counts])
            enrichment=None
            if current_count==0:
                occurrence_rate=-1.0
            else:
                occurrence_rate=codon_chance*overall_count/current_count
            ouput_string=f"{peptide},{occurrence_rate:.4},{current_count}\n"
            if compressing:
                output_file.write(ouput_string.encode())
            else:
               output_file.write(ouput_string)
    output_file.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("uniprot_input_file", help="File containing sequences")
    parser.add_argument("output_filename",
                        help="File to output sequence counts to")
    parser.add_argument("peptide_lengths",
                        help="Length of peptides", type=int)

    # Optional argument flag which defaults to False
    parser.add_argument("-c", "--output_cutoff",
                        help="Dont output sequences apearing more than this many times", action="store", default=False, type=int)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    count_peptides(args.uniprot_input_file, args.output_filename,
                       args.peptide_lengths, args.output_cutoff)