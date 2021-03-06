import os
import sys
import time
import subprocess
from scipy.special import comb
import itertools
import numpy as np
import csv




def netsatsim(Snum, len_p, Len_P):
    input_file = "input"
    mid_file = "mid_file"
    output_file = "output"
    gene_input(Snum, len_p, Len_P)

    with open(input_file, mode='r') as f:
        P = []
        SonP = set()
        color = None
        reader = csv.reader(f)
        for path in reader:
            SonP = SonP.union(path)
            if color is None:
                color = len(path)
            elif color > len(path):
                color = len(path)
            P.append(path)

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


def gene_input(Snum, len_p, Len_P):
    Pset=set()
    Srange=[i for i in range(1, Snum + 1)]
    while(Len_P != len(Pset)):
        Pset.add(tuple(np.random.choice(Srange, len_p, replace=False)))
    with open("input", mode='w') as f:
        writer = csv.writer(f)
        for p in Pset:
            writer.writerow(p)


def gene_mid(P, S, color):
    dimacs = "p cnf {} {}\n".format(color * len(S), int(1 + comb(color, 2)) * len(S) + len(P) * color)
    # one-hot constraint
    for i in range(len(S)):
        line = [i * color + j + 1 for j in range(color)]
        dimacs += ' '.join(map(str, line)) + ' 0\n'
        for k, l in itertools.combinations(line, 2):
            dimacs += "{} {} 0\n".format(-k, -l)

    # output constraint
    Slist=list(S)
    for i in range(len(P)):
        for j in range(color):
            for k in range(len(P[i])):
                dimacs += "{} ".format(int(Slist.index(P[i][k])) * color + j + 1)
            dimacs += '0\n'
    return dimacs



Snummin = 10
Snummax = 30
S_step = 10
len_p = 6
Len_Pmax = 50
iter = 100

Snum=20


with open("result.csv", mode='w') as f:
    for Len_P in range(1, Len_Pmax+1):
        line = []
        for len_p in range(2, 7):
            sum = 0
            for j in range(iter):
                sum+=netsatsim(Snum, len_p, Len_P)
            line.append(sum/iter)
        f.write(str(line).replace("[", "").replace("]", "")+'\n')   #小数点が1桁でしか出力されないためcsvモジュールを使わず
print("simulation is completed!!")
