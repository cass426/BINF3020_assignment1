@dataclass
class Cell:
    score: float
    xp: Optional[int] = None
    yp: Optional[int] = None

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.
    """
    len1 = len(seq1)
    len2 = len(seq2)

    alignment_matrix = [[Cell(0) for x in range(0, len2 + 1)] for y in range(0, len1 + 1)]

    for i in range(1, len1 + 1):
        alignment_matrix[i][0] = Cell(score=i * scoring_function(seq1[i-1], '-'), xp=i-1, yp=0)

    for j in range(1, len2 + 1):
        alignment_matrix[0][j] = Cell(score=j * scoring_function('-', seq2[j-1]), xp=0, yp=j-1)
    
    # doing the scoring bit
    for k in range(1, len1 + 1):
        for l in range(1, len2 + 1):
            diagonal = (alignment_matrix[k-1][l-1].score + scoring_function(seq1[k-1], seq2[l-1]), k-1, l-1)
            up = (alignment_matrix[k-1][l].score + scoring_function(seq1[k-1], '-'), k-1, l)
            left = (alignment_matrix[k][l-1].score + scoring_function('-', seq2[l-1]), k, l-1)

            best = max(diagonal, up, left, key=lambda c: c[0])

            alignment_matrix[k][l] = Cell(score=best[0], xp=best[1], yp=best[2])
    
    # reconstruction
    reconstructed_seq1 = ""
    reconstructed_seq2 = ""
    i = len1
    j = len2
    while (i > 0 or j > 0):
        i_pointer = alignment_matrix[i][j].xp
        j_pointer = alignment_matrix[i][j].yp
        if i_pointer == i:
            reconstructed_seq1 += '-'
            reconstructed_seq2 += seq2[j-1]
            j = j_pointer
            
        elif j_pointer == j:
            reconstructed_seq1 += seq1[i-1]
            reconstructed_seq2 += '-'
            i = i_pointer
        else:
            reconstructed_seq1 += seq1[i-1]
            reconstructed_seq2 += seq2[j-1]
            i = i_pointer
            j = j_pointer

    return reconstructed_seq1[::-1], reconstructed_seq2[::-1], alignment_matrix[len1][len2].score


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    len1 = len(seq1)
    len2 = len(seq2)

    alignment_matrix = [[Cell(0) for x in range(0, len2 + 1)] for y in range(0, len1 + 1)]

    # doing the scoring bit
    for k in range(1, len1 + 1):
        for l in range(1, len2 + 1):
            diagonal = (alignment_matrix[k-1][l-1].score + scoring_function(seq1[k-1], seq2[l-1]), k-1, l-1)
            up = (alignment_matrix[k-1][l].score + scoring_function(seq1[k-1], '-'), k-1, l)
            left = (alignment_matrix[k][l-1].score + scoring_function('-', seq2[l-1]), k, l-1)
            restart = (0, None, None)

            best = max(diagonal, up, left, restart, key=lambda c: c[0])

            alignment_matrix[k][l] = Cell(score=best[0], xp=best[1], yp=best[2])
    
    # reconstruction
    best_score = 0
    best_i, best_j = 0, 0

    if alignment_matrix[k][l].score > best_score:
        best_score = alignment_matrix[k][l].score
        best_i, best_j = k, l
    
    reconstructed_seq1 = ""
    reconstructed_seq2 = ""
    i, j = best_i, best_j

    while alignment_matrix[i][j].xp is not None:
        i_pointer = alignment_matrix[i][j].xp
        j_pointer = alignment_matrix[i][j].yp

        if i_pointer == i:
            reconstructed_seq1 += '-'
            reconstructed_seq2 += seq2[j-1]
        elif j_pointer == j:
            reconstructed_seq1 += seq1[i-1]
            reconstructed_seq2 += '-'
        else:
            reconstructed_seq1 += seq1[i-1]
            reconstructed_seq2 += seq2[j-1]

        i, j = i_pointer, j_pointer

    return reconstructed_seq1[::-1], reconstructed_seq2[::-1], best_score


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

def blosum62_scoring(x, y):
    if x == '-' or y == '-':
        return 8
    return blosum62[x, y]