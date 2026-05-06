# relevant imports
import pytest 
from kmer_analyzer import validate_sequence, update_kmer_count, count_kmers_with_context

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
    assert output["AT"] == 1
    assert output["AT"]["next_chars"] == {"G": 1}
    # result for substring occuring multiple times
    assert output["TG"] == 3
    assert output["TG"]["next_chars"] == {"A": 1, "T": 2}