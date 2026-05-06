# relevant imports
import pytest 
from kmer_analyzer import validate_sequence

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

# test function to validate sequence checking logic
def test_validate_sequence(sequence, k, expected):
    assert validate_sequence(sequence, k) is expected

# basic pytest confirmation
def test_pytest_working():
    assert True
