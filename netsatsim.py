import os
import sys
import subprocess
from scipy.special import comb
import itertools

def gene_mid(P, S, color):
    dimacs = "p cnf {} {}\n".format(color * len(S), int(1 + comb(color, 2)) * len(S) + len(P) * color)
    # one-hot constraint
    for i in range(len(S)):
        line = [i * color + j + 1 for j in range(color)]
        dimacs += ' '.join(map(str, line)) + ' 0\n'
        for k, l in itertools.combinations(line, 2):
            dimacs += "{} {} 0\n".format(-k, -l)

    # output constraint
    for i in range(len(P)):
        for j in range(color):
            for k in range(len(P[i])):
                dimacs += "{} ".format((int(P[i][k].replace('s', '')) - 1) * color + j + 1)
            dimacs += '0\n'
    return dimacs


input_file="sample.txt"
mid_file="mid_file"
output_file="output"
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

while True:
    with open(mid_file, "w") as f:
        f.write(gene_mid(P, S, color))
    subprocess.call(["minisat", mid_file, output_file])
    with open(output_file, "r") as f:
        if f.read=="UNSAT/n":
            color = color - 1
        else :
            break
        if color==0:
            print("could't find satisfiable color sets")
            break