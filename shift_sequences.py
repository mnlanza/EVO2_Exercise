#!/usr/bin/env python3

def read_fasta(filename):
    """Read a FASTA file and return a dictionary of header: sequence pairs."""
    sequences = {}
    current_header = None
    current_sequence = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_sequence)
    
    return sequences

def shift_sequence(sequence, positions):
    """Shift a sequence to the right by n positions, moving start to end."""
    return sequence[-positions:] + sequence[:-positions]

def write_fasta(sequences, filename):
    """Write sequences to a FASTA file."""
    with open(filename, 'w') as f:
        for header, sequence in sequences.items():
            f.write(f"{header}\n")
            # Write sequence in lines of 70 characters
            for i in range(0, len(sequence), 70):
                f.write(f"{sequence[i:i+70]}\n")

def main():
    # Read the original FASTA file
    input_file = "raw_data/gene_variants.fasta"
    sequences = read_fasta(input_file)
    
    # Get the headers in order
    headers = list(sequences.keys())
    
    # Get the 83_S1 sequence that we'll use for shifting
    s1_sequence = sequences['>83_S1']
    
    # Create new sequences dictionary with first two sequences unchanged
    new_sequences = {
        headers[0]: sequences[headers[0]],  # 83_S1
        headers[1]: sequences[headers[1]],  # 83_S2
        headers[2]: shift_sequence(s1_sequence, 1),  # Replace 83_L with shifted 83_S1
        headers[3]: shift_sequence(s1_sequence, 2)   # Replace 83_Stop with shifted 83_S1
    }
    
    # Write the new sequences to a file
    output_file = "shifted_variants.fasta"
    write_fasta(new_sequences, output_file)

if __name__ == "__main__":
    main() 