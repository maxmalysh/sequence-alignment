import numpy
from enum import Enum

class Score:
    match    = +1
    mismatch = -1
    indel    = -3   # gap penalty
    extend   = -1   # gap extension penalty


class Type(Enum):
    Match = 1
    Delete = 2
    Insert = 3


def similarity(first, second):
    return Score.match if first == second else Score.mismatch


def needleman(seq1, seq2, affine=False) -> (str, str):
    m, n = len(seq1), len(seq2)

    # Fill score matrix with initial values
    score = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))
    indel_grid = numpy.zeros(dtype=Type, shape=(m+1, n+1))

    if not affine:
        for i in range(0, m + 1):
            score[i][0] = i * Score.indel
        for j in range(0, n + 1):
            score[0][j] = j * Score.indel

        def score_for(type, i, j):
            score_table = {
                Type.Match  : score[i-1][j-1] + similarity(seq1[i-1], seq2[j-1]),
                Type.Delete : score[i-1][j]   + Score.indel,
                Type.Insert : score[i][j-1]   + Score.indel,
            }
            return score_table[type]

    else:
        score[0][0] = 0
        indel_grid[0][0] = 0

        for i in range(1, m + 1):
            score[i][0] = (i-1) * Score.extend + Score.indel
            indel_grid[i][0] = Type.Delete

        for j in range(1, n + 1):
            score[0][j] = (j-1) * Score.extend + Score.indel
            indel_grid[0][j] = Type.Insert

        def score_for(type, i, j):
            score_table = {
                Type.Match  : score[i-1][j-1] + similarity(seq1[i-1], seq2[j-1]),
                Type.Delete : score[i-1][j]   + Score.extend if indel_grid[i-1][j] == Type.Delete else Score.indel,
                Type.Insert : score[i][j-1]   + Score.extend if indel_grid[i][j-1] == Type.Insert else Score.indel,
            }
            return score_table[type]

    # Calculate scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match  = score_for(Type.Match, i, j)
            delete = score_for(Type.Delete, i, j)
            insert = score_for(Type.Insert, i, j)
            score[i][j] = max(match, delete, insert)

            if affine:
                actions = {
                    match  : Type.Match,
                    delete : Type.Delete,
                    insert : Type.Insert,
                }
                indel_grid[i][j] = actions[score[i][j]]

    # print(score)

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

def smith(seq1, seq2, affine=False) -> (str, str):
    m, n = len(seq1), len(seq2)

    # Generate DP table and traceback path pointer matrix
    score = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))
    pointer = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))

    # initial maximum score in DP table
    max_score = 0

    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i-1][j-1] + similarity(seq1[i-1], seq2[j-1])
            score_up = score[i][j-1] + Score.indel
            score_left = score[i-1][j] + Score.indel
            score[i][j] = max(0,score_left, score_up, score_diagonal)

            # 0 means end of the path
            if score[i][j] == 0:
                pointer[i][j] = 0

            # 1 means trace up
            elif score[i][j] == score_left:
                pointer[i][j] = 1

            #2 means trace left
            elif score[i][j] == score_up:
                pointer[i][j] = 2

            # 3 means trace diagonal
            elif score[i][j] == score_diagonal:
                pointer[i][j] = 3

            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j];

    # initial sequences
    align1, align2 = '', ''

    # indices of path starting point
    i, j = max_i ,max_j

    #traceback, follow pointers
    while pointer[i][j] != 0:
        if pointer[i][j] == 3:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
        elif pointer[i][j] == 1:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1

    return align1[::-1], align2[::-1]

