# Table of Contents

* [What is This?](#what)
* [What Does it Do?](#cool)
* [How Do I Execute?](#execute)
* [File Structure](#structure)
* [Future Goals](#goals)

## What is This? <a name="what"></a>
Implements the Expectation-Maximization (EM) algorithm for transcript quantification from RNA-alignments 
takes in a transcript/alignment file (formatted properly below) and returns estimated 
read-counts for each transcript in the transcriptome.

## What Does it Do? <a name="cool"></a>

Steps:
1. Parse the transcripts into space-efficient data structures
2. Analyze alignment blocks and extract notable information (position, effective length, etc)
3. Calculate relevant probabilities (Pr)
* Pr(Selecting a position along this fragment from which to draw the read)
* Pr(Alignment)
* Pr(Selecting a fragment compatible with the observed length of the one implied by this read)
4. Run the "E" and "M" steps of the EM algorithm until convergence.
5. Observe estimated read-counts for each transcript

## How Do I Execute? <a name="execute"></a>

### (1) Clone
```
git clone https://github.com/RyanDBurke/RNA-seq-quantification-using-the-EM-algorithm.git
```

### (2) Dependencies
* [Python 3.7](https://www.python.org/downloads/)
* [Numpy](http://www.numpy.org/)
* [Numba](https://pypi.org/project/numba/) (Speed-up expensive EM-Steps)
* [SciPy](https://www.scipy.org/)

### (3) Execution
##### then, execute Default or with your own inputs/outputs
```
$ cd data/ 
$ ./squant.py default 
```
```
$ cd data/ 
$ ./squant.py squant --in <input-file>.gz --out <output>
```

### (4) Input File Format
```
number_of_transcripts:int
<transcript_record_1>
<transcript_record_2>
...
<transcript_record_m>
<alignment_block_for_read_1>
<alignment_block_for_read_2>
...
<alignment_block_for_read_n>
```

##### where a transcript record is of the form
```
transcript_name:string <tab> transcript_length:int
```
##### and an alignment block is of the form
```
number_of_alignments_for_read:int
aln1_txp:string <tab> aln1_ori:string <tab> aln1_pos:int <tab> aln1_alignment_prob:double
aln2_txp:string <tab> aln2_ori:string <tab> aln2_pos:int <tab> aln2_alignment_prob:double
...
alnk_txp:string <tab> alnk_ori:string <tab> alnk_pos:int <tab> alnk_alignment_prob:double
```

## File Structure <a name="structure"></a>
    EM-Algorithm
    ├── LICENSE
    ├── README                   
    └── data
        ├── squant.py                       # reads user inputs
        ├── EM.py                           # algorithm in use
        ├── alignments.small.cmsc423.gz     # sample .gz input
        ├── true_counts_small.tsv           # true-counts of sample input
        └── quants.tsv                      # sample output


## Future Goals <a name="goals"></a>
* Improve storage of effective-length computations
* More accurate convergence values to reflect more accurate read-estimations (within 10<sup>-4</sup>)
* Data visualization (spearman correlation plot)
* Allow any input-file. As of now only .gz is accepted

## 
Written for my Bioinformatic Algorithms, Databases, and Tools course
