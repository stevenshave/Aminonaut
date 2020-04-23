class CodonCounter():
    """
    Class which calculates peptide string to codon counts


    AAA is more likely to occcur in randomised codons than WWW. This class
    can be used to calculate the number of times a codon encodes for each
    amino acid in a peptide
    """

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