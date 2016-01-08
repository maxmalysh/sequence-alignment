import argparse
from alignment import NeedlemanAligner, SmithAligner

#   
# Parse input arguments
#
argParser = argparse.ArgumentParser()
argParser.add_argument('filename', type=str, nargs='?',
                       help='Path to the script file.')
argParser.add_argument('--needleman',  dest='needleman',  action='store_true',
                        help='Use Needleman–Wunsch algorithm')
argParser.add_argument('--smith',  dest='smith',  action='store_true',
                        help='Use Smith-Waterman algorithm')
argParser.add_argument('--affine', dest='affine', action='store_true',
                        help='Use affine scoring')

argParser.set_defaults(needleman=False, smith=False)
args = argParser.parse_args()

if args.needleman == args.smith:
    print("You have to choose only one sequence alignment algorithm")
    exit(-1)

if args.filename == None:
    print("You should provide a file with two sequences to be aligned")
    exit(-1)

#
# Read an input file
#
seq = []
with open(args.filename, 'r') as f:
    for line in f:
        line.strip()
        n = -1

        if line.startswith("> seq"):
            n += 1
            seq.append('')
            continue

        if line:
            seq[n] += line.strip().replace(' ', '')

#
# Perform the alignment
#
if args.needleman:
    aligner = NeedlemanAligner(args.affine)
else:
    aligner = SmithAligner(args.affine)

aligned = aligner.align(seq[0], seq[1])

#
# Print output
#
if len(seq[0]) < 20 and len(seq[1]) < 20:
    # Pretty-printing short sequences as in Wikipedia article
    extraws = max(0, len(seq[1])-len(seq[0])), max(0, len(seq[0]) - len(seq[1]))
    extraws = [' ' * (ws+4) for ws in extraws]
    spaces  = ' ' * len(max(extraws[0], extraws[1]))

    print("Sequences%sAligned Sequences" % spaces)
    print("---------%s-----------------" % spaces)
    print("%s %s %s" % (seq[0][0:20], extraws[0], aligned[0][0:20]))
    print("%s %s %s" % (seq[1][0:20], extraws[1], aligned[1][0:20]))
else:
    print(aligned[0], '\n')
    print(aligned[1])



