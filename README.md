# K-mer Analyzer

This project synthesizes pieces of genomic sequences from an input file and determines how often each k-mer exists and the characters that succeed each k-mer across all sequences contained within the input file. This was effectively the final exam for my BIO / DSP 439 class (Big Data Analysis).

## Breakdown

Important features:
  - text file containing example DNA sequences
  - output filename performed through script

Script functionality:
  - validates given sequences
  - counts overlapping k-mers and subsequent sequences
  - combines counts across every sequence in a given file
  - write new lines for each k-mer in a standardized format
      (ex. k = 2 -> TG A:1 T:2)

Running tests:

From the root of the project use this command: 
python kmer_analyzer.py example_sequences.txt 2 output.txt

  - example_sequences.txt: text file with basic example sequences
  - 2: k-mer length
  - output.txt: created for each k-mer summary

For example, the example_sequences.txt may include:
ACGT
ATGTCTGTCTGAA

Then run tests using this command:
python -m pytest

Here is a brief overview of the project structure:
  - DSP439-exam4/                 # parent folder
    - kmer_analyzer.py            # main script
    - example_sequences.txt       # example input file
    - tests/                      # tests folder
      - test_kmer_analyzer.py     # pytest test suite
    - README.md                   # read me file

AI use statement:
I used Perplexity AI to help me troubleshoot bugs and improve the quality of documentation. I also used it to provide test cases for the paremetrize function that I found on the pytest documentation website.
