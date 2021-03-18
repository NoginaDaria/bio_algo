import numpy as np

def create_score_matrix(seq1, seq2, match = 1, mismatch_penalty = -1):
    '''
    create score matrix with custom penalties
    '''
    m = np.full(shape = (len(seq1), len(seq2)), fill_value = mismatch_penalty)
    for i, letter1 in enumerate(seq1):
        for j, letter2 in enumerate(seq2):
            if letter1 == letter2:
                m[i,j] = match
    return m

def fill_scoring_matrix(gap_penalty, scores, n1, n2):
    '''
    fill scoring matrix and store coordinates of parent cells
    '''
    m = np.zeros((n1+1, n2+1))
    parents = {}
    for i in range(1, m.shape[0]):
        for j in range(1, m.shape[1]):
            
            a_b = m[i-1, j-1] + scores[i-1, j-1]
            a_end = m[i-1, j] + gap_penalty
            b_end = m[i, j-1] + gap_penalty

            m[i, j] = np.max([a_b, a_end, b_end, 0])
            
            parent = np.argmax([a_b, a_end, b_end, 0])
            parents[(i, j)] = [(i-1,j-1), (i-1,j), (i,j-1), (0,0)][parent]
    
    return m, parents

def reconstruct_alignment(scoring_matrix, parents, seq1, seq2):
    '''
    uses information from parent cells to build an alignment from scoring matrix
    '''
    max_indeces = np.argwhere(scoring_matrix == np.amax(scoring_matrix))
    
    for max_index in max_indeces:
        curr_index = max_index
        indices = []
        while np.all(curr_index != 0):
            indices.append(curr_index)
            curr_index = parents.get((curr_index[0], curr_index[1]), np.zeros(2))

        seq1_res, seq2_res = '', ''
        indices = indices[::-1]
        for i, ind in enumerate(indices[:-1]):
            seq1_res += seq1[ind[0]] if indices[i][0] != indices[i+1][0] else '-'
            seq2_res += seq2[ind[1]] if indices[i][1] != indices[i+1][1] else '-'
        yield seq1_res, seq2_res, scoring_matrix[max_index[0], max_index[1]]
        
def water(seq1, seq2, match = 1, mismatch = -1, gap_penalty = -1):
    
    '''
    Combines all the computing
    '''
    
    print('Given sequences', seq1, seq2, sep='\n')
    
    scores = create_score_matrix(seq1, seq2, match, mismatch)
    
    scoring_matrix, parents = fill_scoring_matrix(gap_penalty, scores, len(seq1), len(seq2))
    
    alignments = reconstruct_alignment(scoring_matrix, parents, seq1, seq2)
    
    alignments_count = 0
    while True:
        try:
            print('------------')
            align1, align2, score = next(alignments)
            print(f'{align1}\n{align2}\nAlignment score: {score}')
            alignments_count += 1
        except StopIteration:
            print(f'*** Alignments number: {alignments_count} ***')
            break
    return score