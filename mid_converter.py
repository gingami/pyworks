import os
import sys
from scipy.special import comb
input_file="sample.txt"
output_file="mid_file"

if not(os.path.isfile(input_file)):
    print('cannot find "samle.dimacs". ')
    sys.exit()


with open(input_file, mode='r') as f:
    P=[]
    S=set()
    color=None
    for line in f:
        path=line.strip().split(", ")
        S=S.union(path)
        if color is None:
            color=len(path)
        elif color>len(path):
            color=len(path)
        P.append(path)

dimacs="p cnf {} {}".format(color*len(S), )
