import numpy as np

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

    # Observed rates taken by running count_peptides with a peptide length of
    # 1 on UniProt Feb 2018
    aa_to_observed_uniprot_rates={
        '.':1.0,
        'L':0.0965462055734266,
        'A':0.0826589178211278,
        'G':0.0708245377656544,
        'V':0.0687009809306126,
        'E':0.0674303938855079,
        'S':0.0661358493622236,
        'I':0.0592942802510055,
        'K':0.0582300325872045,
        'R':0.0553905654917310,
        'D':0.0546539498083633,
        'T':0.0535108070829884,
        'P':0.0472736671743565,
        'N':0.0405898831521449,
        'Q':0.0393303057739405,
        'F':0.0386161071215318,
        'Y':0.0291769916064555,
        'M':0.0241667625334964,
        'H':0.0227409628455608,
        'C':0.0137794312761490,
        'W':0.0109493679565190,
    }

    def get_codon_occurrence_rate_for_peptide(self, peptide):
        return np.prod([self.aa_to_codon_count[c]/61.0 for c in peptide])
    

    def get_uniprot_observed_occurrence_rate_for_peptide(self, peptide):
        return np.prod([self.aa_to_observed_uniprot_rates[c] for c in peptide])

