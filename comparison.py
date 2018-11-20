import subprocess
from scipy.special import comb
import itertools
import numpy as np
import csv
import random
import functools
import time


def measure(func):
    """
    関数実行時間を計測するラッパー
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start=time.time()
        ret=func(*args, **kwargs)
        return ret, time.time()-start
    return wrapper


def sim(Snum, len_p, Len_P):
    input_file = "input"
    mid_file = "mid_file"
    output_file = "output"
    gene_input(Snum, len_p, Len_P)

    with open(input_file, mode='r') as f:
        P = set()
        reader = csv.reader(f)
        for path in reader:
            P.add(tuple(path))
    return sa_algo(P), netsatsim(P, mid_file, output_file)



def gene_input(Snum, len_p, Len_P):
    Pset=set()
    Srange=[i for i in range(1, Snum + 1)]
    while(Len_P != len(Pset)):
        Pset.add(tuple(np.random.choice(Srange, len_p, replace=False)))
    with open("input", mode='w') as f:
        writer = csv.writer(f)
        for p in Pset:
            writer.writerow(p)



@measure
def sa_algo(Pprev):
    S={Switch(no) for no in genS(Pprev)}
    P=set()
    for pathprev in Pprev:
        path=[]
        for sprev in pathprev:
            for s in S:
                if s.no==sprev:
                    path.append(s)
                    break
        P.add(tuple(path))
    color=0
    colS=set()
    while loop_judge(S, colS, color, P):
        K=gen_U1(P, S,colS, color)
        if C(S, color):
            if len(K)>1:
                K=gen_U2(K, S, color, P)
            if len(K)>1:
                K = gen_U3(K, P)
            if len(K)>1:
                K = gen_U4(K, P)
            if len(K)==0:
                break
        col_s=random.choice(list(K))
        col_s.color=color
        colS.add(col_s)

        if satis_col(P, color, S):
            color+=1

    return color





def gen_U1(P, S,colS, u):
    U1=set()
    max=None
    for s in S.difference(colS):
        length=len(Fa(u, s, P, S))
        if max is None or max<length:
            max=length
    for s in S.difference(colS):
        if len(Fa(u, s, P, S))==max:
            U1.add(s)
    return U1


def gen_U2(K, S, u, P):
    U2=set()
    min = None
    for s in K:
        length = len(Ha(S, u, P, s))
        if min is None or min > length:
            min = length
    for s in K:
        if len(Ha(S, u, P, s)) == min:
            U2.add(s)
    return U2


def gen_U3(K, P):
    U3=set()
    max=None
    for s in K:
        min=None
        for fl in F(P, s):
            length=len(fl)
            if min is None or min > length:
                min=length
        if max is None or max < min:
            max=min
    for s in K:
        min=None
        for fl in F(P, s):
            length=len(fl)
            if min is None or min > length:
                min=length
        if min==max:
            U3.add(s)
    return U3

def gen_U4(K, P):
    U4=set()
    min=None
    for s in K:
        length = len(F(P, s))
        if min is None:
            min = length
        elif min > length:
            min = length
    for s in K:
        if len(F(P, s)) == min:
            U4.add(s)
    return U4







def loop_judge(S, colS, color, P):
    if S==colS:
            return False
    if not color-1<len(list(gen_min_P(P))[0]):
        return False
    return True


def satis_col(P, col, S):
    for path in P:
        if set(path).isdisjoint(C(S, col)):
            return False
    return True



class Switch:
    def __init__(self, no):
        self.no=no
        self.color=None


#フロー集合P上に存在するスイッチ集合を生成
def genS(P):
    S = set()
    for path in P:
        S=S.union(path)
    return S


#スイッチsを経由するフローの集合を生成
def F(P, s):
    F=set()
    for path in P:
        if s in path:
            F.add(path)
    return F


#uで色づけられたスイッチの集合
def C(S, u):
    Su=set()
    for s in S:
        if s.color==u:
            Su.add(s)
    return Su



def Fa(u, s, P, S):
    tmp=s.color
    s.color=u
    Fa=set()
    for path in P:
        if not C(S, u).isdisjoint(set(path)):
            Fa.add(path)
    s.color=tmp
    return Fa


def Ha(S, u, P, s):
    S1=set()
    for swi in C(S, u):
        S1=S1.union(F(P, swi))
    S1=gen_min_P(S1)
    S2=F(P, s)
    return S1.intersection(S2)


@measure
def netsatsim(P, mid_file, output_file):
    color=len(tuple(gen_min_P(P))[0])
    while True:
        with open(mid_file, "w") as f:
            f.write(gen_mid(P, color))
        subprocess.call(["minisat", mid_file, output_file])
        with open(output_file, "r") as f:
            if not(f.read() == "UNSAT\n"):
                return color
            else:
                color = color - 1
        if color == 1:
            print("could't find satisfiable color sets")
            return 1



def gen_min_P(P):
    min=None
    for path in P:
        length=len(path)
        if min is None:
            min = length
        elif min > length:
            min = length
    min_P=set()
    for path in P:
        if len(path)==min:
            min_P.add(path)
    return min_P



def gen_mid(P, color):
    S=genS(P)
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
        sum1 = np.array([0, 0.0])
        sum2 = np.array([0, 0.0])
        for i in range(iter):
            sim1, sim2 = sim(Snum, len_p, Len_P)
            sum1 += np.array(sim1)
            sum2 += np.array(sim2)
        line.extend(sum1 / iter)
        line.extend(sum2 / iter)
        f.write(str(line).replace("[", "").replace("]", "") + '\n')   #小数点が1桁でしか出力されないためcsvモジュールを使わず
print("simulation is completed!!")
