
def levenshtein_distance(seq1, seq2):
    """
    Calculate the levenshtein distance between two nucleotide sequences.
    Used to distinguish between forward and reverse direction subreads.
    """

    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
    if len(seq2) == 0:
        return len(seq1)
    seq1_length = len(seq1) + 1
    seq2_length = len(seq2) + 1
    distance_matrix = [[0] * seq2_length for x in range(seq1_length)]
    for i in range(seq1_length):
        distance_matrix[i][0] = i
        for j in range(seq2_length):
            distance_matrix[0][j]=j
    for i in xrange(1, seq1_length):
        for j in range(1, seq2_length):
            deletion = distance_matrix[i-1][j] + 1
            insertion = distance_matrix[i][j-1] + 1
            substitution = distance_matrix[i-1][j-1]
            if seq1[i-1] != seq2[j-1]:
                substitution += 1
            distance_matrix[i][j] = min(insertion, deletion, substitution)
    return distance_matrix[seq1_length-1][seq2_length-1]