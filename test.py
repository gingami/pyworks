import os
import sys
import subprocess
from scipy.special import comb
import itertools
import numpy as np
import csv

def sim(Snum, len_p, Len_P):
    input_file = "input"
    mid_file = "mid_file"
    output_file = "output"
    gene_input(Snum, len_p, Len_P)

    with open(input_file, mode='r') as f:
        P = set()
        SonP = set()
        reader = csv.reader(f)
        for path in reader:
            path=tuple(path)
            SonP = SonP.union(path)
            P.add(path)
    netsatsim(P, SonP, mid_file, output_file, len(tuple(gen_min_P(P))[0]))



def gene_input(Snum, len_p, Len_P):
    Pset=set()
    Srange=[i for i in range(1, Snum + 1)]
    while(Len_P != len(Pset)):
        Pset.add(tuple(np.random.choice(Srange, len_p, replace=False)))
    with open("input", mode='w') as f:
        writer = csv.writer(f)
        for p in Pset:
            writer.writerow(p)



def gen_min_P(P):
    min=None
    for path in P:
        if min is None:
            min = len(path)
        elif min > len(path):
            min = len(path)
    min_P=set()
    for path in P:
        if len(path)==min:
            min_P.add(path)
    return min_P



def netsatsim(P, SonP, mid_file, output_file, color):
    while True:
        with open(mid_file, "w") as f:
            f.write(gene_mid(P, SonP, color))
        subprocess.call(["minisat", mid_file, output_file])
        with open(output_file, "r") as f:
            if not(f.read() == "UNSAT\n"):
                return color
            else:
                color = color - 1
        if color == 1:
            print("could't find satisfiable color sets")
            return 1



def gene_mid(P, S, color):
    dimacs = "p cnf {} {}\n".format(color * len(S), int(1 + comb(color, 2)) * len(S) + len(P) * color)
    # one-hot constraint
    for i in range(len(S)):
        line = [i * color + j + 1 for j in range(color)]
        dimacs += ' '.join(map(str, line)) + ' 0\n'
        for k, l in itertools.combinations(line, 2):
            dimacs += "{} {} 0\n".format(-k, -l)

    # output constraint
    tup_P=tuple(P)
    Slist=list(S)
    for i in range(len(tup_P)):
        for j in range(color):
            for k in range(len(tup_P[i])):
                dimacs += "{} ".format(int(Slist.index(tup_P[i][k])) * color + j + 1)
            dimacs += '0\n'
    return dimacs




sim(30, 6, 10)

