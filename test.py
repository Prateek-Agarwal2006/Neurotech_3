import sys

def sequence_file(file_path):
    with open(file_path, 'r') as file:
        sequence = file.readline().strip()
    return sequence

def alignment_and_score_calculator(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    m= len(seq1)
    n=len(seq2)
    score_matrix = [[0 for k in range(n+1)] for i in range(m+1)]
    backpropagation_matrix = [[0 for k in range(n+1)] for i in range(m+1)]


    for i in range(1, m+1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(1, n+1):
        score_matrix[0][j] = j * gap_penalty

    
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1][j-1] 
            if seq1[i-1] == seq2[j-1]:
                match += match_score
            else:
                match -= mismatch_penalty

            delete = score_matrix[i-1][j] + gap_penalty

            insert = score_matrix[i][j-1] + gap_penalty

            score_matrix[i][j] = max(match, delete, insert)

            if score_matrix[i][j] == match:
                backpropagation_matrix[i][j] = 'D'
            elif score_matrix[i][j] == delete:
                  backpropagation_matrix[i][j] = 'U'
            else:
                backpropagation_matrix[i][j] = 'L'



    
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        if backpropagation_matrix[i][j] == 'D':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif backpropagation_matrix[i][j] == 'U':
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:  
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, score_matrix[m][n]

def main():
    if len(sys.argv) != 3:
        print("Require python_file.py seq1.txt seq2.txt")
        return

    seq_file1= sys.argv[1]
    seq_file2 = sys.argv[2]
    seq1 = sequence_file(seq_file1)
    seq2 = sequence_file(seq_file2)
    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)

    print("Enter value for match score, mismatch penalty and gap penalty")
    match_score = int(input())
    mismatch_penalty = int(input())
    gap_penalty = int(input())

    alignment1, alignment2, score = alignment_and_score_calculator(seq1, seq2 ,match_score=match_score, mismatch_penalty=mismatch_penalty, gap_penalty=gap_penalty)
    print("Alignment 1:", alignment1)
    print("Alignment 2:", alignment2)
    print("Score:", score)

if __name__ == "__main__":
    main()
