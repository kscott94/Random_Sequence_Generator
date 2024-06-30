# Random_Sequence_Generator
Generate random DNA sequences with a specified GC content

dependencies: Python3, numpy

### Input arguments
NUMBER_SEQUENCES = 3    # how many random sequences do you need?

LENGTH = [200, 300, 150]          # length of sequence(s) in bp. Specify length(s) in python list.

DESIRED_GC_CONTENT = 50 # percent GC content. default=50 percent

ACCEPTABLE_RANGE = 3    # accepatable upper and lower range for GC content. default=3 percent

GC_WINDOW_SIZE = 50     # window size in bp for applying GC content. Window size must be less than the sequence length. Default=50 nt

FORMAT = 'fasta'          # output format for sequence. Options: fasta or tab. default='fasta'

SEED = 321                # seed for reproducibility. default=321
