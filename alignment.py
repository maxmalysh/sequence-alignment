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

class Aligner:
    def align(self, seq1, seq2) -> (str, str):
        raise NotImplementedError

class AffineAligner(Aligner):
    def __init__(self, affine=False):
        self.affine = affine

    def score_generator_for(self, seq1, seq2, scores, indels):
        def score_for(type, i, j):
            if not self.affine:
                score_table = {
                    Type.Match  : scores[i-1][j-1] + similarity(seq1[i-1], seq2[j-1]),
                    Type.Delete : scores[i-1][j]   + Score.indel,
                    Type.Insert : scores[i][j-1]   + Score.indel,
                }
            else:
                score_table = {
                    Type.Match  : scores[i-1][j-1] + similarity(seq1[i-1], seq2[j-1]),
                    Type.Delete : scores[i-1][j]   + Score.extend if indels[i-1][j] == Type.Delete else Score.indel,
                    Type.Insert : scores[i][j-1]   + Score.extend if indels[i][j-1] == Type.Insert else Score.indel,
                }
            return score_table[type]
        return score_for

    def chosen_action(self, match, delete, insert, chosen):
        return {
            match  : Type.Match,
            delete : Type.Delete,
            insert : Type.Insert,
            0      : 0
        }[chosen]


class NeedlemanAligner(AffineAligner):
    def align(self, seq1, seq2) -> (str, str):
        m, n = len(seq1), len(seq2)

        # Fill score matrix with initial values
        score = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))
        indel_grid = numpy.zeros(dtype=Type, shape=(m+1, n+1))

        if not self.affine:
            for i in range(0, m + 1):
                score[i][0] = i * Score.indel
            for j in range(0, n + 1):
                score[0][j] = j * Score.indel
        else:
            score[0][0] = 0
            indel_grid[0][0] = 0

            for i in range(1, m + 1):
                score[i][0] = (i-1) * Score.extend + Score.indel
                indel_grid[i][0] = Type.Delete

            for j in range(1, n + 1):
                score[0][j] = (j-1) * Score.extend + Score.indel
                indel_grid[0][j] = Type.Insert

        score_for = self.score_generator_for(seq1, seq2, score, indel_grid)

        # Calculate scores
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match  = score_for(Type.Match, i, j)
                delete = score_for(Type.Delete, i, j)
                insert = score_for(Type.Insert, i, j)
                score[i][j] = max(match, delete, insert)

                if self.affine:
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

class SmithAligner(AffineAligner):
    def align(self, seq1, seq2) -> (str, str):
        m, n = len(seq1), len(seq2)

        # Generate DP table and traceback path pointer matrix
        score = numpy.zeros(dtype=numpy.int, shape=(m+1, n+1))
        indel_grid = numpy.zeros(dtype=Type, shape=(m+1, n+1))

        score_for = self.score_generator_for(seq1, seq2, score, indel_grid)

        # initial maximum score
        max_score = 0

        # Calculate scores
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match  = score_for(Type.Match, i, j)
                delete = score_for(Type.Delete, i, j)
                insert = score_for(Type.Insert, i, j)

                # 0 means end of the path
                score[i][j] = max(0, match, delete, insert)
                if self.affine:
                    indel_grid[i][j] = {
                        match  : Type.Match,
                        delete : Type.Delete,
                        insert : Type.Insert,
                        0      : 0
                    }[score[i][j]]

                if score[i][j] >= max_score:
                    max_i = i
                    max_j = j
                    max_score = score[i][j]

        # print(score)

        align1, align2 = '', ''

        # indices of path starting point
        i, j = max_i, max_j

        # Traceback
        while score[i][j] != 0:
            if score[i][j] == score_for(Type.Match, i, j):
                align1 += seq1[i-1]
                align2 += seq2[j-1]
                i -= 1
                j -= 1
            elif score[i][j] == score_for(Type.Insert):
                align1 += '-'
                align2 += seq2[j-1]
                j -= 1
            elif score[i][j] == score_for(Type.Delete):
                align1 += seq1[i-1]
                align2 += '-'
                i -= 1

        while i > 0:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1

        while j > 0:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

        # print("Score is ", score[max_i][max_j])

        return align1[::-1], align2[::-1]

