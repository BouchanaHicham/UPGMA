# ----------------------------------- [ Profile Matrix ] ------------------------------------

def Profile_Matrix(sequences):
    if not sequences:
        print("Sequences Are Empty Fam")

    # Check that all sequences have the same length
    seq_length = len(sequences[0])
    if not all(len(seq) == seq_length for seq in sequences):
        print("Sequences Do Not Have The Same Length so they are not aligned")

    # Initialize profile matrix with zeros
    profile_matrix = {nucleotide: [0] * seq_length for nucleotide in "ACGT-"}

    # Update profile matrix based on aligned sequences
    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        for nucleotide in "ACGT-":
            frequency = column.count(nucleotide) / len(column)  # 1/3 for example
            profile_matrix[nucleotide][i] = round(frequency, 2)  # Round to 2 decimal place

    return profile_matrix

# ----------------------------------- [ Check If Tansition ] ------------------------------------
def is_transition(nucleotide1, nucleotide2):
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    return (nucleotide1 in purines and nucleotide2 in purines) or (nucleotide1 in pyrimidines and nucleotide2 in pyrimidines)
    # This will return True if Transition

# ------------------------------------ [ Score Matrix ] ------------------------------------

def calculate_score(nucleotide_j, position):
    MNuc = 2
    MGap = 1
    MMs = 1
    MMv = -1
    gap_penalty = -3
    score = 0

    for nucleotide, frequencies in profile_matrix.items():
        if nucleotide == nucleotide_j:
            # Set the weight based on the nucleotide type
            if nucleotide_j == '-':
                weight = MGap
            else:
                weight = MNuc
        elif nucleotide == '-' or nucleotide_j == '-':
            weight = gap_penalty
        elif is_transition(nucleotide, nucleotide_j):
            weight = MMs
        else:
            weight = MMv
        # Calculate the score for the current nucleotide and position
        score += frequencies[position] * weight

    return score

# ------------------------------------ [ Align Sequences ] ------------------------------------


def align_sequences(seqs, T):
    m, n = len(T), len(seqs[0])
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    direction_matrix = [[''] * (n + 1) for _ in range(m + 1)]

    gap_penalty = -3
    # Initialize the first column of score_matrix and direction_matrix
    for i in range(1, m + 1):
        score_matrix[i][0] = round(score_matrix[i - 1][0] + gap_penalty, 3)
        direction_matrix[i][0] = '↑'

    # Initialize the first row of score_matrix and direction_matrix
    for j in range(1, n + 1):
        score_matrix[0][j] = round(score_matrix[0][j - 1] + gap_penalty, 3)
        direction_matrix[0][j] = '←'

    for i in range(1, m + 1):
        for j in range(1, n + 1):

            delete = round(score_matrix[i - 1][j] + calculate_score('-', j - 1), 3)
            insert = round(score_matrix[i][j - 1] + calculate_score('-', j - 1), 3)

            match = round(score_matrix[i - 1][j - 1] + calculate_score(T[i - 1], j - 1), 3)

            max_score = max(match, delete, insert)

            score_matrix[i][j] = max_score
            if max_score == match:
                direction_matrix[i][j] += '↖'
            if max_score == delete:
                direction_matrix[i][j] += '↑'
            if max_score == insert:
                direction_matrix[i][j] += '←'

    print("Score Matrix:")
    for row in score_matrix:
        print([round(val, 3) for val in row])
    print("\nDirection Matrix:")
    for row in direction_matrix:
        print(row)

    all_alignments = []

    def backtrack(i, j, seq_align):
        nonlocal all_alignments

        if i == 0 and j == 0:
            all_alignments.append(seq_align[::-1])
            return

        if '↖' in direction_matrix[i][j]:
            backtrack(i - 1, j - 1, seq_align + T[i - 1])
        if '↑' in direction_matrix[i][j]:
            backtrack(i - 1, j, seq_align + '-')
        if '←' in direction_matrix[i][j]:
            backtrack(i, j - 1, seq_align + '-')

    backtrack(m, n, '')

    return all_alignments, score_matrix[m][n]



aligned_sequences = ["AC-GT", "AC-GT", "GCCAT"]
T = "ACG"
# ------------ [ Profile Matrix ] ------------

print("Profile Matrix:")
profile_matrix = Profile_Matrix(aligned_sequences)
# Print the profile matrix
for nucleotide, frequencies in profile_matrix.items():
    print(f"{nucleotide}: {frequencies}")

# ------------ [ Align Sequences ] ------------
aligned_result, s = align_sequences(aligned_sequences, T)
print("Aligned Result:")
for r in aligned_sequences:
    print(r)
print(aligned_result[0])
print("Alignment Score:", s)

# Double check
