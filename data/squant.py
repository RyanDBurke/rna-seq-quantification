#! /usr/bin/env python3

"""
* Our main driver for our program.
* This will read our command-line inputs
* execute with ./squant.py squant --in alignments.cs423 --out quants.tsv --eqc
*               [0]        [1]    [2]  [3]              [4]   [5]        [6]
* --eqc is optional
"""

import sys
import EM

def squant():

    # default execution
    if (sys.argv[1].lower() == "default"):
        input = "alignments_small.cs423.gz"
        output = "quants.tsv"

    else:

        # determine if a command was entered
        if len(sys.argv) < 6:
            print("please enter a valid tool")
        else:

            squant = sys.argv[1]
            in_ = sys.argv[2]
            input = sys.argv[3]
            out_ = sys.argv[4]
            output = sys.argv[5]

            # check to see if flags were correctly typed
            if (not squant == "squant") or (not in_ == "--in") or (not out_ == "--out"):
                print('\033[31m' + "one of your flags was mistyped.")
                print('\033[39m')
                return

    import time
    start = time.time()

    # call EM
    EM.EM_Algorithm(input, output)

    # print time
    time = (time.time()-start)/60
    print('runtime ~ ', round(time, 1), 'minutes')

    # Nice little print-out
    print('\033[31m' + "You can find results in " + output)

    # reset to default color
    print('\033[39m') 

"run our program"
squant()

