from witch_msa.helpers.alignment_tools import Alignment
from pyhmmer import easel

# helper function to convert an alignment object to TextMSA object
def alignmentToTextMSA(aln, name):
    sequences = []
    for taxon, seq in aln.items():
        sequences.append(easel.TextSequence(name=taxon.encode(), sequence=seq))
    msa = easel.TextMSA(name=name.encode(), sequences=sequences)
    return msa

# helper function to obtain alphabet given molecule type
def moleculeToAlphabet(molecule):
    alphabet = None
    if molecule == 'amino':
        alphabet = easel.Alphabet.amino()
    elif molecule == 'dna':
        alphabet = easel.Alphabet.dna()
    elif molecule == 'rna':
        alphabet = easel.Alphabet.rna()
    else:
        raise ValueError(f"alphabet {molecule} is not amino, dna, or rna")
    return alphabet
