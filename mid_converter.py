import os
import sys
from scipy.special import comb
import itertools

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

dimacs="p cnf {} {}\n".format(color*len(S), (int(1+comb(len(S), 2))+len(P))*color)
#one-hot constraint
for i in range(color):
    line=[i+1+j*color for j in range(len(S))]
    dimacs+=' '.join(map(str, line))+'\n'
    for k, l in itertools.combinations(line, 2):
        dimacs+="{} {}\n".format(-k, -l)

#output constraint
for i in range(len(P)):
    for j in range(color):
        for k in range(len(P[i])):
            dimacs+="{} ".format((int(P[i][k].replace('s', ''))-1)*color+j+1)
        dimacs+='\n'
print(dimacs)