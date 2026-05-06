# relevant imports
import pytest 
import sys
from kmer_analyzer import (
    validate_sequence,
    update_kmer_count,
    count_kmers_with_context,
    write_results_to_file,
    main,
)

# basic pytest confirmation
def test_pytest_working():
    assert True

# dynamic test functionality through parametrize
@pytest.mark.parametrize(
    "sequence,k,expected",
    [
        ("ACGT", 2, True),
        ("ATGTCTGTCTGAA", 2, True),
        ("A", 2, False),
        ("", 1, False),
        ("AC1T", 2, False),
        ("1234", 2, False),
    ],
)

#############################
## validate sequence tests ##
#############################

# testing basic genomic 
def test_validate_sequence(sequence, k, expected):
    assert validate_sequence(sequence, k) is expected

##############################
## update k-mer count tests ##
##############################

# confirm equence is sufficient length and only holds valid characters
def test_update_correct_kmer_counts():
    # empty dictionary to hold kmer data
    kmer_data = {}
    # recorded observations of k-mer and following character
    output = update_kmer_count(kmer_data, "AC", "G")
    # total count for record 
    assert output == {
        "AC": {
            "count": 1,
            "next_chars": {"G": 1}
        }
    }

# confirm sequence follows constraints with same next character
def test_update_kmer_counts_with_existing_count():
    # previously observed dictionary to hold kmer data
    kmer_data = {
        "AC": {
            "count": 1,
            "next_chars": {"G": 1}
        }
    }
    # record same structured observation
    output = update_kmer_count(kmer_data, "AC", "G")
    # increased total count of recorded observation
    assert output == {
        "AC": {
            "count": 2,
            "next_chars": {"G": 2}
        }
    }

# confirm sequence follows constraints with new character
def test_update_kmer_counts_with_new_next_char():
    # same dictionary to hold kmer data with existing record
    kmer_data = {
        "AC": {
            "count": 1,
            "next_chars": {"G": 1}
        }
    }
    # record different next character for same k-mer
    output = update_kmer_count(kmer_data, "AC", "T")
    # increased total count of recorded observation
    assert output == {
        "AC": {
            "count": 2,
            "next_chars": {"G": 1, "T": 1}
        }
    }

###############################
## k-mers with context tests ##
###############################

# confirm correct processing with sequencing lacking sufficient length
def test_sequence_lacking_context():
    # sequence test example shrorter than k + 1
    sequence = "A"
    k = 2
    # record k-mers and following character counts
    output = count_kmers_with_context(sequence, k)
    # resulting empty dictionary due to insufficient sequence length
    assert output == {}

# confirm correct processing with sequencing of sufficient length
def test_sequence_with_sufficient_context():
    # sequence with one k-mer and following character
    sequence = "ACGT"
    k = 2
    # record k-mers and following character counts
    output = count_kmers_with_context(sequence, k)
    # resulting substring and following character counts
    assert output == {
        "AC": {
            "count": 1,
            "next_chars": {"G": 1}
        },
        "CG": {
            "count": 1,
            "next_chars": {"T": 1}
        }
    }

# confirm correct processing with k-mers and following characters observed multiple times
def test_sequence_with_repeated_kmers():
    # sequence from assignment example
    sequence = "ATGTCTGTCTGAA"
    k = 2
    # record k-mers and following character counts
    output = count_kmers_with_context(sequence, k)
    # result for substring occuring one time
    assert output["AT"]["count"] == 1
    assert output["AT"]["next_chars"] == {"G": 1}
    # result for substring occuring multiple times
    assert output["TG"]["count"] == 3
    assert output["TG"]["next_chars"] == {"T": 2, "A": 1}

###################################
## written results to file tests ##
###################################

# test to create sample file with expected test format
def test_write_results_to_file(tmp_path):
    # example k-mer data to write to file
    kmer_data = {
        "AC": {
            "count": 1,
            "next_chars": {"G": 1}
        },
        "TG": {
            "count": 3,
            "next_chars": {"T": 2, "A": 1}
        }
    }
    # create temporary file path for output
    output_file = tmp_path / "output.txt"
    # call to test writing results to file
    write_results_to_file(kmer_data, output_file)
    # read file contents
    contents = output_file.read_text()
    # expected file contents with sorted k-mers and next characters
    assert contents == "AC G:1\nTG A:1 T:2\n"

#########################
## main function tests ##
#########################

def test_main_counts_acros_multiple_sequences(tmp_path):
    # create temporary input file with two part sequences
    input_file = tmp_path / "input.txt"
    input_file.write_text("ACGT\nATGTCTGTCTGAA\n")
    # create temporary output file path
    output_file = tmp_path / "output.txt"
    # save original argument to hold after test
    original_argv = sys.argv
    # simulate script arguments
    sys.argv = ["kmer_analyzer.py", str(input_file), "2", str(output_file)]
    # run actual main function to process input file and write results
    try:
        main()
    finally:
        sys.argv = original_argv
    # read file contents
    contents = output_file.read_text()
    # expected file contents with sorted k-mers and next characters from both sequences
    assert contents == ("AC G:1\nAT G:1\nCG T:1\nCT G:2\nGT C:2\nTG A:1 T:2\n")