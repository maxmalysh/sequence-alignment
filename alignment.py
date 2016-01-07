import numpy
from enum import Enum

class Score:
    match    = +1
    mismatch = -1
    indel    = -1   # gap penalty
    extend   =  0   # gap extension penalty


class Type(Enum):
    Match = 1
    Delete = 2
    Insert = 3


def similarity(first, second):
    return Score.match if first == second else Score.mismatch


def needleman(seq1, seq2, affine=False):
    m, n = len(seq1), len(seq2)

    # Fill score matrix with initial values
    score = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))

    for i in range(0, m + 1):
        score[i][0] = i * Score.indel
    for j in range(0, n + 1):
        score[0][j] = j * Score.indel

    def score_for(type, i, j):
        if type is Type.Match:
            return score[i-1][j-1] + similarity(seq1[i-1], seq2[j-1])
        elif type is Type.Delete:
            return score[i-1][j] + Score.indel
        elif type is Type.Insert:
            return score[i][j-1] + Score.indel

    # Calculate scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match  = score_for(Type.Match, i, j)
            delete = score_for(Type.Delete, i, j)
            insert = score_for(Type.Insert, i, j)
            score[i][j] = max(match, delete, insert)

    print(score)

    # Traceback and compute the alignment
    # starting from the bottom right cell
    align1, align2 = '', ''
    i,j = m,n

    while i > 0 and j > 0: # end touching the top or the left edge
        current = score[i][j]

        # A diagonal arrow represents a match or mismatch,
        # either way this means each letter corresponds to another letter.
        if current == score_for(Type.Match, i, j):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1

        # If there is a horizontal arrow there will be two columns for one row in the alignment
        # and hence a gap in the side string. This gap is after the letter in the row.
        elif current == score_for(Type.Delete, i, j):
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1

        # If there is a vertical arrow there will be two rows for one column in the alignment
        # and hence a gap in the top string. This gap is after the letter in the column
        elif current == score_for(Type.Insert, i, j):
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

        #print("{0}\n{1}\n".format(align1, align2))

    # If sequences have different length, we have to complete the longest one,
    # so we move further to the topmost left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1

    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    return align1[::-1], align2[::-1]


def smith(seq1, seq2, affine=False):
    return seq1
