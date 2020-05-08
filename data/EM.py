"""
* EM-Algorithm
"""

from numba import jit
import gzip
import scipy.stats
import numpy as np
import os

"used to speed up effective length calculations"
def memoizer(f):
    memolist = [-1] * 1000
    def _f(l):
            if memolist[l] == -1:
                    memolist[l] = f(l)
            return memolist[l]
    return _f

"effective length for shorter transcripts"
@memoizer
def get_eff_len_short(l):
    mu = 200
    sd = 25
    d = scipy.stats.norm(mu, sd)
    # the (discrete) distribution up to 200
    p = np.array([d.pdf(i) for i in range(l+1)])
    # re-normalize so this is a proper distribution
    p /= p.sum()
    # the expected value of a distribution f is \sum_{i=0}^{max_i} f(i) * i
    cond_mean = np.sum([i * p[i] for i in range(len(p))])
    return l - cond_mean

"returns effective length"
def effective_length(l):
    mu = 200
    return l - mu if l >= 1000 else get_eff_len_short(l)


"returns (p2 * p3 * p4)"
def calculateP(l, position, eff, origin, prob):

    # true length
    length = l

    # pos
    pos = position

    # calc p2
    p2 = 1 / eff

    # calc p3
    p3 = 0
    d = scipy.stats.norm(200, 25)
    D = d.cdf
    if origin == "f":
        p3 = D(length - pos)
    else:
        p3 = D(pos + 100)

    # calc p4
    p4 = prob

    # returns p2 * p3 * p4
    return (p2 * p3 * p4)

"Runs our EM Algorithm Steps"
@jit(nopython=True)
def runEM(numTranscripts, t_index, p_list, block, est_reads, it):

    # make n 1/M for all transcripts
    n = np.full(numTranscripts, (1 / numTranscripts))

    # next estimate for n in the next iteration of the EM
    n_next = np.zeros(numTranscripts)

    # number of alignment blocks
    each_block = len(block) - 1

    "E-Step"
    while True:

        # have all transcripts converged?
        converged = True
        it = it + 1

        # for each alignment block
        for b in range(each_block):

            # normalized denomiator
            normal = 0.0

            # for each record in this current block
            for record in range(block[b], block[b + 1]):

                # index of current transcript
                index_of_transcript = t_index[record]

                # p1*p2*p3*p4 of current transcript
                p = n[index_of_transcript] * p_list[record]

                # update normal
                normal = normal + p
                    
            # for each record in this current block
            for record in range(block[b], block[b + 1]):

                # index of current transcript
                index_of_transcript = t_index[record]

                # p1*p2*p3*p4 of current transcript
                p = n[index_of_transcript] * p_list[record]

                # update estimated num of reads for each transcript
                est_reads[index_of_transcript] = est_reads[index_of_transcript] + p
                
        "M-Step"
        # calculate total est num of reads
        N = 0
        for estimate in est_reads:
            N = N + estimate

        # update n
        index_of_transcript = 0
        for estimate in est_reads:

            # new n for this transcript
            newN = estimate / float(N)

            # update n_next
            n_next[index_of_transcript] = newN

            # increment
            index_of_transcript = index_of_transcript + 1

        # have all transcripts converged?
        for c in range(len(n)):

            # compare n and n'
            old = n[c]
            new = n_next[c]

            # avoid divding by zero
            if new > 0:

                # check if its within 10^(-2) or 0.01%
                check = abs(new - old) / new
                if not check <= 0.01:
                    converged = False
                    break

        # if it converged then let's break
        if converged:
            return est_reads, it

        # otherwise, let's go back to the E-Step
        else:
            n = n_next
            n_next = np.zeros(numTranscripts)

"main() of EM Algorithm"
def EM_Algorithm(in_, out_):
    txt = False

    print("\nhang tight, this might take a while...")

    # determine type of file extension
    fileType = os.path.splitext(in_)[1]

    # open based on whether its .txt or .gz
    if fileType == ".gz":

        # open .gz input file
        in_ = gzip.open(in_, "r")

    else:

        # open regular file
        in_ = open(in_, "r")

        if fileType == ".txt":
            txt = True

    # output file
    out_ = out_

    # each line of input file
    lines = in_.readlines()

    # map of transcripts
    transcripts = dict()

    # reverse index of transcripts, helps access transcript by index
    transcriptList = []

    # number of transcript (M + 1)
    if txt:
        numTranscripts = int(lines[0]) + 1
    else:
        numTranscripts = int(lines[0].decode('utf-8')) + 1

    # index of transcript
    index = 0

    # input transcripts into map
    for i in range(1, numTranscripts):

        # get each transcript
        if txt:
            transcript = lines[i]
        else:
            transcript = lines[i].decode('utf-8')

        transcript = transcript.split()

        # get their name, length, and effective length
        name = transcript[0]
        length = int(transcript[1])
        effLength = effective_length(length)

        # place into our map
        transcripts[name] = [length, index, effLength]

        # place into our reverse index list
        transcriptList.append(name)

        # increment index
        index += 1

    print("(1/5) transcript parsing done")
    print("(2/5) alignment blocks done")

    # we want to start on the line after all the transcripts
    line_number = numTranscripts

    # End-Of-File
    EOF = len(lines)

    # holds indexes of each transcript
    t_index = []

    # (p2*p3*p4) for each transcript
    p_list = []

    # length of each alignment block
    block = [0]

    # read till the EOF
    while (line_number < EOF):

        # num of alignments in this read
        if txt:
            num_of_align = int((lines[line_number].split())[0])
        else:
            num_of_align = int((lines[line_number].decode('utf-8').split())[0])

        # add alignments to read
        for a in range(line_number + 1, line_number + num_of_align + 1):

            # split the alignment record
            if txt:
                align = lines[a].split()
            else:
                align = lines[a].decode('utf-8').split()

            # name
            name = align[0]

            # origin
            origin = align[1]

            # pos
            pos = int(align[2])

            # prob
            p4 = float(align[3])

            # length, index, eff
            l, index, eff = transcripts[name]

            # add index to t_index
            t_index.append(index)

            # calculate p2*p3*p4
            p = calculateP(l, pos, eff, origin, p4)

            # (p2*p3*p4)
            p_list.append(p)

        # update block length
        block.append(len(p_list))    

        # increment to next alignment block
        line_number += num_of_align + 1

    if txt:
        numTranscripts = int(lines[0])
    else:
        numTranscripts = int(lines[0].decode('utf-8'))
    
    # cast to numpy arrays and creating EM parameters
    t_index = np.array(t_index)
    p_list = np.array(p_list)
    block = np.array(block)
    est_reads = np.zeros(numTranscripts)
    it = 0

    print("(3/5) running EM")
    final_read_est, it = runEM(numTranscripts, t_index, p_list, block, est_reads, it)
    print("(4/5) EM finsihed")
    s = ""
    s += "name\teff-length\test_frags\n"

    print("(5/5) writing to quants.tsv \n")
    it = '{:,}'.format(it)
    print(it + " iterations of EM ran \n")
    i = 0
    for transcript in transcriptList:

        # serialize
        name = transcript
        length = (transcripts[transcript])[2]
        est = final_read_est[i]

        # to string
        s += name + "\t"
        s += str(round(float(length), 0)) + "\t"
        s += str(est) + "\n"

        # increment
        i += 1


    # write to our output
    f = open(out_, "w")
    f.write(s)
    f.close()





        





                




