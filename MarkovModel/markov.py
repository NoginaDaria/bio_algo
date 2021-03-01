from collections import defaultdict
from tqdm import tqdm
import numpy as np

# contains : sequences, Bernoully, Markov1st, Markov2nd, Markov3rd

# !!!the module can read only from fasta file with one sequence and header!!!

def sequences(file, n = 1000, l = 4000):
    '''
    Generator creating n sequences of a given length l from fasta file
    '''
    with open(file) as f:
        f.readline()
        curr_seq = ''
        n_count = 0
        while n_count < n:
            seq = f.readline()
            if not seq:
                break
            else:
                seq = seq.strip()
            if l - len(curr_seq) > len(seq):
                curr_seq += seq
            else:
                curr_seq += seq[:l - len(curr_seq)]
            if len(curr_seq) == l:
                n_count += 1
                yield curr_seq
                curr_seq = ''
                
class Bernoully():
    '''
    Bernoully model
    '''
    def __init__(self):
        self.order = 0
        
    def fit(self, file, l = 4000, n = 1000):
        self.l = l
        seq_gen = sequences(file, self.l, n)
        self.transition_matrix = defaultdict(lambda: 0)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for i in s[self.order:]:
                self.transition_matrix[i] += 1
        
        dict_size = sum(list(self.transition_matrix.values()))
        for amino_acid in self.transition_matrix:
            self.transition_matrix[amino_acid] = self.transition_matrix[amino_acid] / dict_size
            
    def gen(self, l=None):
        if not l:
            l = self.l
        res = ''.join(np.random.choice(list(self.transition_matrix.keys()), l, list(self.transition_matrix.values())))
        
        return res
    
                
class Markov1st():
    '''
    First-order Markov model
    '''
    def __init__(self):
        self.order = 1
                
    def fit(self, file, l = 4000, n = 1000):
        self.l = l
          
        # create transition matrix of co-occurence of words
        
        self.transition_matrix = defaultdict(lambda: defaultdict(lambda: 0))
        
        # also create first-occurence dict
        self.first_letters = {0: defaultdict(lambda: 0)}

        seq_gen = sequences(file, self.l, n)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            self.first_letters[0][s[0]] += 1
            for i in range(1, len(s)):
                self.transition_matrix[s[i-1]][s[i]] += 1
                    
        # normalize transition matrix
        
        # create nucleotide dict with their occurence
        
        seq_gen = sequences(file, self.l, n)
        amino_dict = defaultdict(lambda: 0)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for i in s[self.order:]:
                amino_dict[i] += 1
                
        # normalizing
        
        for word in tqdm(self.transition_matrix):
            for tr_word in self.transition_matrix[word]:
                num_occ = sum(self.transition_matrix[word].values())
                self.transition_matrix[word][tr_word] = self.transition_matrix[word][tr_word] / num_occ
                
        # normalize first-occurence dict
        first_sum = sum(self.first_letters[0].values())
        for letter in self.first_letters[0]:
            self.first_letters[0][letter] = self.first_letters[0][letter] / first_sum
        
    
    def gen(self, l=None):
        if not l:
            l = self.l
        
        res = np.random.choice(list(self.first_letters[0].keys()), 1, p=list(self.first_letters[0].values()))[0]
        
        for i in range(self.order, l):
            d = self.transition_matrix[res[i-1]]
            p = list(d.values())
            res += np.random.choice(list(d.keys()), 1, p=p)[0]
            
        return res
    
    
    
class Markov2nd():
    '''
    Second-order Markov model
    '''
    def __init__(self):
        self.order = 2
                
    def fit(self, file, l = 4000, n = 1000):
        self.l = l 
        
        # create transition matrix of co-occurence of words
        # first word in the dict goes first
        
        self.transition_matrix = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
        
        # also create first-occurence dict
        
        self.first_letters = {0: defaultdict(lambda: 0),
                              1: defaultdict(lambda: 0)}
        
        seq_gen = sequences(file, self.l, n)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for order in range(self.order):
                self.first_letters[order][s[order]] += 1
        
            for i in range(self.order, len(s)):
                self.transition_matrix[s[i-2]][s[i-1]][s[i]] += 1
                    
        # normalize transition matrix
        
        # create nucleotide dict with their occurence
        seq_gen = sequences(file, self.l, n)
        amino_dict = defaultdict(lambda: 0)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for i in s[self.order:]:
                amino_dict[i] += 1
                
        for order in range(self.order):
            first_sum = sum(self.first_letters[order].values())
            for letter in self.first_letters[order]:
                self.first_letters[order][letter] = self.first_letters[order][letter] / first_sum
                
        # normalizing
        
        for word in tqdm(self.transition_matrix):
            for tr_word_1 in self.transition_matrix[word]:
                for tr_word_2 in self.transition_matrix[word][tr_word_1]:
                    num_occ = sum(self.transition_matrix[word][tr_word_1].values())
                    self.transition_matrix[word][tr_word_1][tr_word_2] = self.transition_matrix[word][tr_word_1][tr_word_2] / num_occ
    
    def gen(self, l=None):
        if not l:
            l = self.l
        
        res = ''
        for order in range(self.order):
            p = np.array(list(self.first_letters[order].values()))
            res += np.random.choice(list(self.first_letters[order].keys()), 1, p=p/p.sum())[0]
        
        first, second = res[0], res[1]
        for i in range(self.order, l):
            d = self.transition_matrix[first][second]
            p = list(d.values())
            res += np.random.choice(list(d.keys()), 1, p=p)[0]
            first, second = second, res[-1]
            
        return res
    
    
    
class Markov3rd():
    '''
    Third-order Markov model
    '''
    def __init__(self):
        self.order = 3
                
    def fit(self, file, l = 4000, n = 1000):
        self.l = l
        
        # create transition matrix of co-occurence of words
        # first word in the dict goes first
        
        self.transition_matrix = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0))))
        
        # also create first-occurence dict
        
        self.first_letters = {0: defaultdict(lambda: 0),
                              1: defaultdict(lambda: 0),
                              2: defaultdict(lambda: 0)}
        
        seq_gen = sequences(file, self.l, n)
        
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for order in range(self.order):
                self.first_letters[order][s[order]] += 1
        
            for i in range(self.order, len(s)):
                self.transition_matrix[s[i-3]][s[i-2]][s[i-1]][s[i]] += 1
                    
        # normalize transition matrix

        # create nucleotide dict with their occurence
        seq_gen = sequences(file, self.l, n)
        
        amino_dict = defaultdict(lambda: 0)
        for _ in tqdm(range(n)):
            s = next(seq_gen)
            for i in s[self.order:]:
                amino_dict[i] += 1
                
        for order in range(self.order):
            first_sum = sum(self.first_letters[order].values())
            for letter in self.first_letters[order]:
                self.first_letters[order][letter] = self.first_letters[order][letter] / first_sum
                
        # normalizing
        
        for word in tqdm(self.transition_matrix):
            for tr_word_1 in self.transition_matrix[word]:
                for tr_word_2 in tr_word_1:
                    for tr_word_3 in tr_word_2:
                        num_occ = sum(self.transition_matrix[word][tr_word_1][tr_word_2].values())
                        self.transition_matrix[word][tr_word_1][tr_word_2][tr_word_3] = self.transition_matrix[word][tr_word_1][tr_word_2][tr_word_3] / num_occ
    
    def gen(self, l=None):
        if not l:
            l = self.l
        
        res = ''
        for order in range(self.order):
            p = np.array(list(self.first_letters[order].values()))
            res += np.random.choice(list(self.first_letters[order].keys()), 1, p=p/p.sum())[0]
        
        first, second, third = res[0], res[1], res[2]
        for i in range(self.order, l):
            d = self.transition_matrix[first][second][third]
            p = list(d.values())
            res += np.random.choice(list(d.keys()), 1, p=p)[0]
            first, second, third = second, res[-2], res[-1]
            
        return res