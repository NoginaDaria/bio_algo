import numpy as np
import markov

def alphabet_size(file):
    '''
    Get size of thex alphabet
    Hardcore because of the aim to generalize the approach
    '''
    alphabet = set()
    with open(file) as f:
        f.readline()
        lines = [i.strip() for i in f.readlines()]
        for line in lines:
            alphabet.update(set([letter for letter in line]))
    return len(alphabet)

def get_param_num(model, file):
    '''
    Number of independent parameters of a given model
    '''
    k = model.order
    A = alphabet_size(file)
    
    res = 0
    for i in range(k):
        res += A ** (k + 1)
        
    return res

def AIC(k, L):
    '''
    Akaike Information Criterion 
    Log of likelihood was precalculated
    '''
    return k * 2 - 2 * L

def BIC(k, L, n):
    '''
    Bayesian Information Criterion
    Log of likelihood was precalculated
    '''
    return k * np.log(n) - 2 * L

def calc_model_score(model, file, n, l):
    '''
    Calculate AIC and BIC of the given model
    '''
    k = get_param_num(model, file)
    
    seq_gen = markov.sequences(file=file, n=n, l=l)
    
    # calculate Likelihood
    
    Ls = []
    
    for _ in range(n):
        L = 0
        seq = next(seq_gen)
    
        # process first letters probability
        
        for order in range(model.order):
            first_nums = model.first_letters[order]
            p = first_nums.get(seq[order])
            assert p <= 1, 'Wrong probability'
            L += np.log(first_nums.get(seq[order]) / sum(first_nums.values()))
           
        # other letters probability
        for letter in seq[model.order:]:
            prob = model.transition_matrix[letter]
            for order in range(model.order):
                prob = prob.get(letter, 'not found')
                assert prob != 'not found', 'Letter combination not found'
                    
            L += np.log(prob)
        Ls.append(L)

    return AIC(k, sum(Ls)), BIC(k, sum(Ls), n)