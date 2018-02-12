import numpy as np
import scipy
import itertools as it
from sys import argv

file_1 = open("mixed2.xyz", "r")


def printMat(a):
    rows = a.shape[0]
    cols = a.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            print("% 2.5f  " % a[i, j]),
        print
    print


class atomlistclass(object):
    def __init__(self):
        self.num = 0
        self.sym = ''
        self.mass = 0
        self.name = ''
        self.radius = ''

    def read_atomlist(self, num, sym, mass, name, r):
        self.num = float(num)
        self.sym = sym
        self.mass = float(mass)
        self.name = name
        if r == 'NA':
            self.radius = r
        else:
            self.radius = float(r)


atomfile = open("atom.csv", "r")
atomlist = {}
for i, line in enumerate(atomfile):
    if i > 0:
        (at_num, at_sym, at_mass, at_name, r) = line.split()
        atomlist[i] = atomlistclass()
        atomlist[i].read_atomlist(at_num, at_sym, at_mass, at_name, r)


def massfromatom(findit):
    for i in atomlist:
        if atomlist[i].sym == findit:
            mass = atomlist[i].mass
    return mass


def covradfromatom(findit):
    for i in atomlist:
        if atomlist[i].sym == findit:
            r = atomlist[i].radius
    return r


class atom(object):
    def __init__(self):
        self.index = 0
        self.mol = ''
        self.x = 0
        self.y = 0
        self.z = 0

    def read_atom(self, i, n, x, y, z, mass):
        self.index = i
        self.mol = n
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass

    def cov_rad(self, r):
        self.radius = r


A = {}
for i, line in enumerate(file_1):
    if i > 1:
        (n, x, y, z) = line.split()
        X = float(x)
        Y = float(y)
        Z = float(z)
        mass = massfromatom(n)
        r = covradfromatom(n)
        A[i - 2] = atom()
        A[i - 2].read_atom(i, n, X, Y, Z, mass)
        A[i - 2].cov_rad(r)


def bond_length(A):
    n = len(A)
    D = np.zeros((n, n))
    Adj = np.zeros((n, n))
    Adj2 = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            dist = ((A[i].x - A[j].x)**2 + (A[i].y - A[j].y)**2 +
                    (A[i].z - A[j].z)**2)**0.5
            D[i, j] = dist

            if i != j and dist < A[i].radius + A[j].radius:
                Adj[i, j] = 1
                Adj2[i, j] = (A[i].radius + A[j].radius) / (
                    A[i].mass * A[j].mass)**0.5

    return D, Adj, Adj2


D, Adj, Adj2 = bond_length(A)
printMat(Adj)

P = np.dot(Adj, Adj)
printMat(P)

print np.diag(P)
printMat(Adj)
printMat(Adj2)

P2 = np.dot(Adj2, Adj2)
printMat(P2)
for i in range(0, len(A)):
    print A[i].mol, " ", P[i, i], " ", P2[i, i]
