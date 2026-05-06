import sys

def validate_sequence(sequence, k):
    if len(sequence) < k:
        return False
    for nucleotide in sequence:
        if nucleotide in '1234567890':
            return False
    return True

def update_kmer_count(kmer_data, kmer, next_char):
    if kmer not in kmer_data:
        # small revision
        # initiatize k_mer as if it is not seen prior to observation
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    kmer_data[kmer]['count'] += 1
    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    kmer_data = {}
    
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            f.write(f"{kmer} {next_char_str}\n")


def main():
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")

    # small addition
    # dictionary holding cumulative k-mer data counts
    all_kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            ## kmer_data = count_kmers_with_context(sequence, k) ## not included now
            
        # large addition
        # specifically count current k-mer seuqneces
        sequence_kmer_data = count_kmers_with_context(sequence, k)
        # merge current sequence counts into cumulative dictionary
        for kmer, data in sequence_kmer_data.items():
            if kmer not in all_kmer_data:
                all_kmer_data[kmer] = {'count': 0, 'next_chars': {}}
            # add counts from current sequence to cumulative dictionary
            all_kmer_data[kmer]['count'] += data['count']
            # add next character counts from current sequence to cumulative dictionary
            for next_char, freq in data['next_chars'].items():
                if next_char not in all_kmer_data[kmer]['next_chars']:
                    all_kmer_data[kmer]['next_chars'][next_char] = 0
                all_kmer_data[kmer]['next_chars'][next_char] += freq
        # small revision
        # output progress after each sequence is processed
        write_results_to_file

if __name__ == '__main__':
    main()