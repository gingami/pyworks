import numpy as np
import csv
Len_S=50
len_p=10
Len_P=100
with open("input", mode='w') as f:
    writer=csv.writer(f)
    for i in range(Len_P):
        p=list(np.random.choice([i for i in range(1, Len_S + 1)], len_p, replace=False))
        writer.writerow(p)