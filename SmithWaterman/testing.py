from random import randint, choice
def random_seq(n=15):
    return ''.join([choice('ATCG') for _ in range(n)])

def change_seq(seq, n=3, mutation=True, ends=False):
    mutated_pos = []
    while len(mutated_pos) < n + 1:
        
        ends_pos = choice([randint(0, len(seq)//5), randint(len(seq)-len(seq)//5, len(seq)-1)])
        pos = randint(1, len(seq)-1) if not ends else ends_pos
        
        if mutation:
            if pos not in mutated_pos:
                seq = seq[:pos] + choice('ATCG') + seq[pos+1:]
                print(len(seq))
        else:
            seq = seq[:pos] + seq[pos+1:]
            
        mutated_pos.append(pos)
    return seq