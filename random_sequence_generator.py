import numpy as np

### set your parameters, then hit run #############################################################

NUMBER_SEQUENCES = 3    # how many random sequences do you need?
LENGTH = [200, 300, 150]          # length of sequence(s) in bp. Specify length(s) in python list.
DESIRED_GC_CONTENT = 50 # percent GC content.
ACCEPTABLE_RANGE = 3    # accepatable upper and lower range for GC content. 
GC_WINDOW_SIZE = 50     # window size in bp for applying GC content.
FORMAT = 'fasta'          # output format for sequence. Options: fasta or tab
SEED = 321                # seed for reproducibility.

##################################################################################################

# Input parameters to string (will print to console)
params = f'NUMBER_SEQUENCES = {NUMBER_SEQUENCES}   # how many random sequences do you need?\n\
LENGTH = {LENGTH}           # length of sequence(s). Specify length(s) in python list\n\
DESIRED_GC_CONTENT = {DESIRED_GC_CONTENT} # percent GC content\n\
ACCEPTABLE_RANGE = {ACCEPTABLE_RANGE}    # accepatable upper and lower range for GC content\n\
GC_WINDOW_SIZE = {GC_WINDOW_SIZE}     # window size in bp for applying GC content\n\
FORMAT = {FORMAT}          # output format for sequence (fasta or tab)\n\
SEED = {SEED}                # seed for reproducibility\n'


# set the seed
np.random.seed(SEED)

# helpers
def get_random_sequence(length):
    """
    Generate random sequence of specified length.
    :param length: int, desired sequence length
    :return: chr, random sequence
    """
    alphabet = ['A', 'T', 'C', 'G']
    output = ''
    for i in range(length):
        pos = np.random.randint(len(alphabet))
        output = output + alphabet[pos]
    return(output)


def get_GC(seq):
    """
    calculate GC content of specified sequence.
    :param seq: chr, sequence
    :return: int, GC content percent
    """
    GCs = 0
    for i in seq:
        if i =="G" or i == "C":
            GCs += 1
    percent_GC_content = round(100*(GCs / len(seq)))
    return(percent_GC_content)


def get_GC_window(seq, window_size = 50):
    """
    Calculate the highest and lowest GC content in specified window size
    :param seq: chr, entire sequence
    :param window_size: int, nucleotide length of window size
    :return: int, percent GC content of the lowerst and highest GC content windows.
    """
    windows = len(seq) - window_size
    position1 = 0
    position2 = window_size
    highest_windowed_GC = 0
    lowest_windowed_GC = 100

    for i in range(windows):
        windowed_GC_itteration = get_GC(seq[position1:position2])

        if windowed_GC_itteration > highest_windowed_GC:
            highest_windowed_GC = windowed_GC_itteration
        if windowed_GC_itteration < lowest_windowed_GC:
            lowest_windowed_GC = windowed_GC_itteration

        position1 += 1
        position2 += 1

    return(highest_windowed_GC, lowest_windowed_GC)


def get_good_window(length, desired_GC, acceptable_range):
    """
    When encountering a window with percent GC content outside of desired range, randomly regenerate the sequence
    within that windown until the GC content falls within desired range.
    :param length: integer, windown length
    :param desired_GC: integer, percent GC content
    :param acceptable_range: int, percent GC content flexibility
    :return: chr, new sequence with desired GC content
    """
    seq = get_random_sequence(length)
    GC_content = get_GC(seq)

    while GC_content > (desired_GC + acceptable_range) or \
            GC_content < (desired_GC - acceptable_range):
        seq = get_random_sequence(length)
        GC_content = get_GC(seq)

    return(seq)


def main(length, desired_GC=50, acceptable_range=2, window_size=50):
    """
    Generate random nucleotide sequence(s) with a specific percent GC content
    :param length: list, list of ints defning the nucleotide length of each sequence
    :param desired_GC: int, desired percent GC content. default=50 %
    :param acceptable_range: int, how far higher or lower the GC content may be. default=2 %
    :param window_size: How many nucleotides within a sliding window to check GC content. default=50 nucleotides
    :return: chr, fasta or tab formated nucleotide sequences.
    """
    rand_sequence = get_random_sequence(length)
    position1 = 0
    position2 = window_size
    itter_position = len(rand_sequence)

    while position2 <= itter_position:
        current_window = rand_sequence[position1:position2]
        CG_current_window = get_GC(current_window)

        current_window_check = False
        if (desired_GC - acceptable_range) <= CG_current_window <= (desired_GC + acceptable_range):
            current_window_check = True
            position1 += 1
            position2 += 1

        if current_window_check == False:
            window_replacement = get_good_window(window_size, desired_GC, acceptable_range)
            rand_sequence = rand_sequence[:position1] + window_replacement + rand_sequence[position2:]
            position1 = 0
            position2 = window_size


    final_length = len(rand_sequence)
    final_GC = get_GC(rand_sequence)
    final_GC_window = get_GC_window(rand_sequence,window_size)
    GC_windowed_h = final_GC_window[0]
    GC_windowed_l = final_GC_window[1]

    return(rand_sequence, final_length, final_GC, GC_windowed_h, GC_windowed_l)


# get random sequence(s) and print output to console
if __name__ == "__main__":
    if FORMAT != 'fasta' and FORMAT != 'tab':
        print("Error: Please specify output format. Options: FORMAT = 'fasta' or 'tab'")
        exit(0)

    if type(LENGTH) != list:
        print("LENGTHs should be specified in a python list encompased in brackets and separated by commas.")
        print("e.g. LENGTH = [250,195,300] ")
        print("e.g. LENGTH = [250] ")
        exit(0)

    if len(LENGTH) != NUMBER_SEQUENCES:
        print(f'you listed {len(LENGTH)} lengths in LENGTH but indicated {NUMBER_SEQUENCES} random sequences to be gernated.   ')
        print("Please specify one length for every sequence.")
        exit(0)

    print(params)  # print input parameters
    print_header = True  # for tab output, print column headers only once

    for i in range(NUMBER_SEQUENCES):
        generate_random_sequence = main(LENGTH[i],DESIRED_GC_CONTENT,ACCEPTABLE_RANGE,GC_WINDOW_SIZE)
        random_sequence = generate_random_sequence[0]
        final_length = generate_random_sequence[1]
        GC_random_sequence = generate_random_sequence[2]
        highest_windowed_GC = generate_random_sequence[3]
        lowest_windowed_GC = generate_random_sequence[4]

        # print output
        if FORMAT == 'fasta':
            header = f'>sequence{i+1} length:{final_length}, GC content:{GC_random_sequence}%, actual windowed GC content range:{lowest_windowed_GC}-{highest_windowed_GC}%'
            print(header)
            print(random_sequence)
        elif FORMAT == 'tab':
            if print_header == True:
                print(f'sequence number\trandom sequence\tlength\tGC content\tactual windowed GC content range')
                print_header = False
            print(f'{i+1}\t{random_sequence}\t{final_length}\t{GC_random_sequence}%\t{lowest_windowed_GC}-{highest_windowed_GC}%')
