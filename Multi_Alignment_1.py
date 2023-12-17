import numpy as np

def build_profile_matrix(sequences):
    symbols = set(''.join(sequences))
    profile_matrix = {symbol: [0] * len(sequences[0]) for symbol in symbols}

    for seq in sequences:
        for i, symbol in enumerate(seq):
            profile_matrix[symbol][i] += 1
    print("profile_matrix",profile_matrix)
    return profile_matrix

def align_sequences_with_profile(template, profile_matrix, sequences):
    aligned_sequences = []

    for seq in sequences:
        alignment = []

        for i in range(len(template)):
            nucleotide_template = template[i]
            nucleotide_seq = seq[i]

            if profile_matrix[nucleotide_template][i] >= profile_matrix[nucleotide_seq][i]:
                alignment.append(nucleotide_template)
            else:
                alignment.append(nucleotide_seq)

        aligned_sequences.append(''.join(alignment))

    return aligned_sequences

def dynamic_programming_alignment(template, profile_matrix, sequences, match_score=2, gap_penalty=-3):
    N, M = len(sequences), len(template)

    score_matrix = np.zeros((N + 1, M + 1), dtype=int)
    direction_matrix = np.zeros((N + 1, M + 1), dtype=int)

    for i in range(1, N + 1):
        score_matrix[i, 0] = score_matrix[i - 1, 0] + gap_penalty

    for j in range(1, M + 1):
        score_matrix[0, j] = score_matrix[0, j - 1] + gap_penalty

    for i in range(1, N + 1):
        for j in range(1, M + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if template[j - 1] == sequences[i - 1][j - 1] else 0)
            gap_up = score_matrix[i - 1, j] + gap_penalty
            gap_left = score_matrix[i, j - 1] + gap_penalty

            score_matrix[i, j] = max(match, gap_up, gap_left)

            if score_matrix[i, j] == match:
                direction_matrix[i, j] = 1  # Diagonal
            elif score_matrix[i, j] == gap_up:
                direction_matrix[i, j] = 2  # Up
            else:
                direction_matrix[i, j] = 3  # Left

    aligned_sequences = []
    i, j = N, M

    while i > 0 or j > 0:
        if direction_matrix[i, j] == 1:  # Diagonal
            aligned_sequences.append(sequences[i - 1])
            i -= 1
            j -= 1
        elif direction_matrix[i, j] == 2:  # Up
            aligned_sequences.append('-' * M)
            i -= 1
        else:  # Left
            aligned_sequences.append(template)
            j -= 1

    aligned_sequences.reverse()
    print("direction_matrix",direction_matrix)
    return aligned_sequences

# Example usage:
template_sequence = 'ACG'
aligned_sequences_input = ['AC-GT', 'AC-GT', 'GCCAT']

profile_matrix = build_profile_matrix(aligned_sequences_input)
aligned_sequences_output = align_sequences_with_profile(template_sequence, profile_matrix, aligned_sequences_input)

print("Aligned Sequences with Profile:")
for seq in aligned_sequences_output:
    print(seq)

dynamic_programming_aligned_sequences = dynamic_programming_alignment(template_sequence, profile_matrix, aligned_sequences_input)

print("\nAligned Sequences using Dynamic Programming:")
for seq in dynamic_programming_aligned_sequences:
    print(seq)
