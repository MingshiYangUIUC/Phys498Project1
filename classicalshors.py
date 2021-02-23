#imports
import numpy as np
from sympy import isprime
from math import gcd
import time
import scipy
import scipy.sparse
from numpy.random import normal,uniform,randint
from matplotlib import pyplot as plt
import sys
from fractions import Fraction

def checkeasy(N):
    easy = 0
    isEven = 0
    isPrime = 0
    isXpa = 0
    if N < 3:
        easy = 1
    if easy == 0 and N % 2 == 0:
        easy = 1
        isEven = 1
    if easy == 0 and isprime(N) == True:
        easy = 1
        isPrime = 1
    if easy == 0:
        log2N = int(np.log10(N)/np.log10(2))
        for i in range(2,log2N):
            x = N**(1/i)
            if np.abs(x - np.round(x,0))<1e-9:
                easy = 1
                isXpa = 1
            if easy == 1:
                break
    return easy,isEven,isPrime,isXpa



def xN_gcd(N):
    x = randint(2,np.sqrt(N)+1)
    while gcd(x,N) != 1:
        x = randint(2,np.sqrt(N))
    return x

def xr1modN(x,N):
    r = 1
    while ((x**r)%N) != (1 % N):
        r += 1
    return r
#print(xN_gcd(122211))

#print(xr1modN(52,122211))

def factor(N):
    if checkeasy(N)[0] == 1:
        return 'easy'
    else:
        ans = []
        passed = 0
        ntry = 0
        while passed == 0 and ntry < 10:
            r = 1
            while (r % 2) != 0:
                x = xN_gcd(N)
                r = xr1modN(x,N)
                ntry += 1
            sol = [gcd(int((x**int(r/2)-1)%N),N),gcd(int((x**int(r/2)+1)%N),N)]
            for elem in sol:
                if elem != 1 and elem != N:
                    passed = 1
                    if elem not in ans:
                        ans.append(elem)
        return ans

def factorv2(N):
    if checkeasy(N)[0] == 1:
        return 'easy'
    else:
        ans = []
        passed = 0
        ntry = 0
        while passed == 0 and ntry < 100:
            r = 1
            while (r % 2) != 0 or r == 0:
                x = xN_gcd(N)
                #r = xr1modN(x,N)
                phases = np.array(np.linalg.eig(xNUnitary(x,N))[0])
                eigenvals = np.log(phases)/2/np.pi/(1.j)
                value = np.abs(eigenvals[np.random.randint(0,len(eigenvals))])
                print(value)
                r = Fraction(value).limit_denominator(100).denominator
                print(r)
                ntry += 1
            sol = [gcd(int((x**int(r/2)-1)%N),N),gcd(int((x**int(r/2)+1)%N),N)]
            for elem in sol:
                if elem != 1 and elem != N:
                    passed = 1
                    if elem not in ans:
                        ans.append(elem)
        for d in ans:
            f = N/d
            if f not in ans:
                ans.append(int(f))
        if ans == []:
            print('!!!!!!!')
            print(N,x,r,sol,value)
            print('-------')
        return ans


def nbits(N):
    return len('{0:b}'.format(N))



def xNUnitary(x,N):
    n = int(np.ceil((np.log10(N)/np.log10(2)))) #size of qubits, means n wires?
    #print(n)
    matrix = np.identity(2**n)*0
    for i in range(2**n):
        if i < N:
            a = i
            b = (i*x)%N
            #print(a,b)
            matrix[b][a] += 1
            #matrix[i][i] *= ((i*x)%N)
        else:
            matrix[i][i] += 1
    return matrix


#print(xN_gcd(10))
#print(xr1modN(3,10))
#print(xNUnitary(3,10))
#print(np.linspace(0,15,16) @ xNUnitary(3,10))

#phases = np.array(np.linalg.eig(xNUnitary(3,10))[0])

#eigenvals = np.log(phases)/2/np.pi/(1.j)
#value = np.abs(eigenvals[np.random.randint(0,len(eigenvals))])
#print(value)
#print(eigenvals*xr1modN(3,10))
#print((Fraction(0.6).limit_denominator(100).denominator))
'''
for i in range(100):
    v = factorv2(i)
    if type(v) != str:
        print(i,v)
        pass
'''
print(factorv2(21))
'''
print(xN_gcd(23))
print(xr1modN(3,10))
basis = '0'
while len(basis) < 5:
    basis = '0' + basis
print(basis)'''