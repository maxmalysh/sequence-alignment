import argparse

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

argParser.set_defaults(needleman=False, smith=False)
args = argParser.parse_args()

if args.needleman == args.smith:
    print("You have to choose only one sequence alignment algorithm")
    exit(-1)

if args.filename == None:
    print("You should provide a file with two sequences to be aligned")
    exit(-1)

