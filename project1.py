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

#time parameters
T0 = 0
T1 = 0
T2 = 0
T3 = 0

#classical shor's algorithm functions
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
    ntry = 0
    while gcd(x,N) != 1 and ntry < 100:
        x = randint(2,np.sqrt(N)+1)
        ntry += 1
    #if ntry == 30:
    #    print(f'faila{N}')
    return x
def xr1modN(x,N):
    r = 1
    ntry = 0
    while ((x**r)%N) != (1 % N) and ntry < 100:
        r += 1
        ntry += 1
    #if ntry == 30:
    #    print(f'failb{N}')
    return r
def factor(N): #this one has been modified to accept "Easy numbers" for a part of question
    if checkeasy(N)[0] == 2:
        return 'easy'
    else:
        ans = []
        passed = 0
        ntry = 0
        while passed == 0 and ntry < 100:
            r = 1
            while (r % 2) != 0 and ntry < 100:
                x = xN_gcd(N)
                r = xr1modN(x,N)
                ntry += 1
            sol = [gcd(int((x**int(r/2)-1)%N),N),gcd(int((x**int(r/2)+1)%N),N)]
            for elem in sol:
                if elem != 1 and elem != N and ntry<100:
                    passed = 1
                    if elem not in ans:
                        ans.append(elem)
        if ntry == 100 and passed == 0:
            print(f'failc{N}')
        if (1 not in getsol(x,r,N)) and (N not in getsol(x,r,N)):
            return ans,x,r
        else:return(f'failed{N}',ans,x,r,getsol(x,r,N))
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
                u = xNUnitary(x,N)
                lbd = np.array(np.linalg.eig(u)[0])
                phases = np.log(lbd)/2/np.pi/(1.j)
                el = np.abs(phases)
                minim = 1
                
                for k in el:
                    if k < minim and k > 0.0001:
                        minim = k
                    else:
                        pass
                
                #print(minim)
                rand = el[np.random.randint(0,len(el))]
                value = np.abs(rand)
                r = Fraction(value).limit_denominator(100).denominator
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
def xNUnitary(x,N):
    n = int(np.ceil((np.log10(N)/np.log10(2)))) #size of qubits, means n wires?
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
def getsol(x,r,N):
    sol = [gcd(int((x**int(r/2)-1)%N),N),gcd(int((x**int(r/2)+1)%N),N)]
    return sol

#basic functions
def IntToBinary(integer):
    return int('{0:b}'.format(integer))
def IntToBinary_OLD(integer): 
    if integer == 0:
        return 0
    else:
        save=integer
        digits=int(np.ceil(np.log10(integer)/np.log10(2)))
        numb=0
        if save == 2**digits:
            numb = 10**digits
        else:
            for i in range(digits):
                val = int(np.floor(integer / (2**(digits-i-1))))
                integer = int(integer % (2**(digits-i-1)))
                numb += int(val*(10**(digits-i-1)))
        return numb
def BinaryToInt(binary): # input is string of a number or integer
    return int(str(binary),2)
def BinaryToInt_OLD(binary): # input is string of a number or integer
    binstr=str(binary)
    numb=0
    j=len(binstr)
    for i in range(j):
        numb += int(binstr[j-i-1])*(2**(i))
    return numb
def PrettyPrintBinary(state,printline=0):
    string = '( '
    for element in state:
        if printline == 0:
            string += (str(element[0])+' |'+element[1]+'> + ')
        else:
            string += ('\n'+str(element[0])+' |'+element[1]+'> + ')
    newvaling= string[:-3]
    newvaling += ' )'
    print(newvaling)
    pass
def PrettyPrintInteger(state,printline=0):
    string = '( '
    for element in state:
        if printline == 0:
            string += (str(element[0])+' |'+str(BinaryToInt(element[1]))+'> + ')
        else:
            string += ('\n'+str(element[0])+' |'+str(BinaryToInt(element[1]))+'> + ')
    newvaling= string[:-3]
    newvaling += ' )'
    print(newvaling)
    pass
def StateToVec(state):
    values = []
    indexes = []
    for element in state:
        values.append(element[0])
        indexes.append(BinaryToInt(element[1]))
    srted = list(indexes)
    srted.sort()
    maxindex = (srted[-1])
    maxindex += 1
    nwires = len(state[0][1])
    vector = np.zeros(2**nwires,dtype=complex)
    i=0
    for index in indexes:
        k = values[i]
        vector[index] = k
        i+=1
    return vector
def VecToState(vector):
    state = []
    lenstr = np.log10(len(vector))/np.log10(2)
    for i in range(len(vector)):
        value = vector[i]
        index = str(IntToBinary(i))
        if value != 0:
            while len(index) < lenstr:
                index = '0' + index
            state.append((value,index))
    return state

#Matrixes: Hadamard, phase, CNOT. identity matrix: np.identity(size)
def Hadamard():
    H = np.ndarray((2,2),dtype=complex)*0+1
    H[1][1] = -1
    return H/(np.sqrt(2))
def Phase(theta): #theta is in radians
    P = np.ndarray((2,2),dtype=complex)*0
    P[0][0] = 1
    P[1][1] = np.exp(theta*1.j)
    return P
def CNOT(ctrl,other):
    C = np.identity(4,dtype=complex)
    if ctrl < other:
        C[2][2] = 0
        C[2][3] = 1
        C[3][2] = 1
        C[3][3] = 0
    elif ctrl > other:
        C[1][1] = 0
        C[1][3] = 1
        C[3][1] = 1
        C[3][3] = 0
    return C

#tensor product function (list of matrices)
def tensorMe(listOfMatrices):
    matri = listOfMatrices[0]
    for i in range(len(listOfMatrices)):
        if i != 0:
            matri = np.kron(matri,listOfMatrices[i])
    return matri
#arrays with parameters to build the unitary matrix

def HadamardArray(i, k):#i_th wire on totally k wires, start from 0...
    # this should apply Hadamard to wire i out of k wires
    size=2**(k+1)
    #myMatrix=np.zeros((size,size),dtype=complex)
    matrilist = []
    for j in range(k+1):#when j==i tensorproduct with hadamard otherwise by I
        if i != j:
            matrilist.append(np.identity(2,dtype=complex))
        else:
            matrilist.append(Hadamard())
    myMatrix = tensorMe(matrilist)
    return myMatrix
def PhaseArray(i,theta,k,):#i_th wire on totally k wires, start from 0...
    # this should apply Hadamard to wire i out of k wires
    size=2**(k+1)
    matrilist = []
    for j in range(k+1):#when j==i tensorproduct with phase otherwise by I
        if i != j:
            matrilist.append(np.identity(2,dtype=complex))
        else:
            matrilist.append(Phase(theta))
    myMatrix = tensorMe(matrilist)
    return myMatrix
def CNOTArray(ctrl,other,k):#control, other wire on totally k wires, start from 0...
    firstwire = min(ctrl,other)
    matrilist = []
    j = 0
    while j < (k+1):#when j == smaller of the two tensorproduct with cnot, and i += 1, otherwise by I
        if j != firstwire:
            matrilist.append(np.identity(2,dtype=complex))
            j += 1
        else:
            matrilist.append(CNOT(ctrl,other))
            j += 2
    myMatrix = tensorMe(matrilist)
    return myMatrix

#read file and generate unitary matrix or gate list. decompose of gate is available in the read_decompose
def ReadCircuit(fileName):
    myInput_lines=open(fileName).readlines()
    myInput=[]
    numberOfWires=int(myInput_lines[0])
    for line in myInput_lines[1:]:
        myInput.append(line.split())
    return (numberOfWires,myInput)
def Read_decompose(fileName):
    myInput_lines=open(fileName).readlines()
    myInput=[]
    numberOfWires=int(myInput_lines[0])
    inputVec = []
    for line in myInput_lines[1:]:
        gate = line.split()
        #print(gate)
        if gate[0] == 'INITSTATE':
            if gate[1] == 'FILE':
                inputVec = ReadVec(gate[2])
            if gate[1] == 'BASIS':
                inputVec = StateToVec([[(1+0.j),gate[2][1:-1]]])
        elif gate[0] == 'H' or gate[0] == 'P' or gate[0] == 'CNOT':
            myInput.append(gate)
        elif gate[0] == 'FUNC' or gate[0] == 'CFUNC':
            myInput.append(gate)
        elif gate[0] == 'NOT' or gate[0] == 'RZ' or gate[0] == 'CRZ' or gate[0] == 'CPHASE' or gate[0] == 'SWAP':
            decomp = []
            inwire = gate[1]
            if gate[0] == 'NOT':
                decomp.append(['H',inwire])
                decomp.append(['P',inwire,'np.pi'])
                decomp.append(['H',inwire])
            elif gate[0] == 'RZ':
                decomp.append(['P',inwire,eval(gate[2])/2])
                decomp.append(['H',inwire])
                decomp.append(['P',inwire,'np.pi'])
                decomp.append(['H',inwire])
                decomp.append(['P',inwire,eval(gate[2])/(-2)])
                decomp.append(['H',inwire])
                decomp.append(['P',inwire,'np.pi'])
                decomp.append(['H',inwire])
            elif gate[0] == 'CRZ':
                otherwire = gate[2]
                decomp.append(['P',otherwire,eval(gate[3])/2])
                decomp.append(['CNOT',inwire,otherwire])
                decomp.append(['P',otherwire,eval(gate[3])/(-2)])
                decomp.append(['CNOT',inwire,otherwire])
            elif gate[0] == 'CPHASE':
                otherwire = gate[2]
                decomp.append(['P',inwire,eval(gate[3])/2])
                decomp.append(['P',otherwire,eval(gate[3])/2])
                decomp.append(['CNOT',inwire,otherwire])
                decomp.append(['P',otherwire,eval(gate[3])/-2])
                decomp.append(['CNOT',inwire,otherwire])
            elif gate[0] == 'SWAP':
                otherwire = gate[2]
                decomp.append(['CNOT',otherwire,inwire])
                decomp.append(['CNOT',inwire,otherwire])
                decomp.append(['CNOT',otherwire,inwire])
            else:
                pass
            for elem in decomp:
                myInput.append(elem)
        else:
            myInput.append(gate)
        #print(myInput)
    return (numberOfWires,myInput,inputVec)
def UnitaryMatrix(nWires,myInput):#return the matrix, and 0 or 1 for whether MEASURE
    RevGateList = []
    matrices = []
    measure = 0
    if myInput[-1] == ['MEASURE']:
        measure = 1
        del myInput[-1]
    for i in range(len(myInput)-1,-1,-1):
        RevGateList.append(myInput[i])
    for gate in RevGateList:
        if gate[0] == 'H':
            matrices.append(HadamardArray(int(gate[1]),nWires-1))
        elif gate[0] == 'P':
            matrices.append(PhaseArray(int(gate[1]),float(gate[2]),nWires-1))
        elif gate[0] == 'CNOT':
            matrices.append(CNOTArray(int(gate[1]),float(gate[2]),nWires-1))
    Unitary = matrices[0]
    j = 0
    for matrix in matrices:
        if j != 0:
            Unitary = np.dot(Unitary,matrix)
        j += 1
    return Unitary,measure
def UnitaryFromFile(filename):#return the matrix, and 0 or 1 for whether measure as like (Matrix,0)
    InputData = Read_decompose(filename)
    nWires = InputData[0]
    myInput = InputData[1]
    inputVec = InputData[2]
    return UnitaryMatrix(nWires,myInput),inputVec,nWires
def ReadVec(fileName):
    State_lines=open(f'project1/{fileName}').readlines()
    vecs = np.zeros(len(State_lines),dtype=complex)
    for i in range(len(State_lines)):
        nums = State_lines[i].split()
        vecs[i] = complex(nums[0]) + complex(nums[1]+'j')
    return vecs

#measure function to obtain probabilities and output
def MeasureProb(state):
    binaryandprob=[]
    for i in range(len(state)):
        data = state[i]
        if data[0] != 0:
            binaryandprob.append([data[1],(np.absolute(data[0]))**2])
    return binaryandprob
def ProbToOut(probs):
    p = np.random.uniform()
    out = 'NA'
    for i in range(len(probs)):
        if i == 0:
            pbase = 0
        else: 
            pbase = pcap
        pcap = (pbase + probs[i][1])
        if pbase < p and p <= pcap:
            out = probs[i][0]
        #print(out,probs[i],p,pbase,pcap)
    return out


#begin simulator.1.b, different operation method!
def UnitaryOperation(nWires,myInput,vector):#The operation!, and retain 0 or 1 for whether MEASUREprint(myInput)
    GateList = []
    matrices = []
    measure = 0
    if myInput[-1] == ['MEASURE']:
        measure = 1
        del myInput[-1]
    for i in range(len(myInput)):
        GateList.append(myInput[i])
    for gate in GateList:
        if gate[0] == 'H':
            matrices.append(HadamardArray(int(gate[1]),nWires-1))
        elif gate[0] == 'P':
            matrices.append(PhaseArray(int(gate[1]),float(gate[2]),nWires-1))
        elif gate[0] == 'CNOT':
            matrices.append(CNOTArray(int(gate[1]),int(gate[2]),nWires-1))
    for matrix in matrices:
        vector = np.dot(matrix,vector)
    return vector,measure

#Begin simulator 1c functions needed to be modified to sparsed matrix
#tensor product function (list of matrices)
def tensorMeSP(listOfMatrices):
    matri = listOfMatrices[0]
    for i in range(len(listOfMatrices)):
        if i != 0:
            matri =scipy.sparse.kron(matri,listOfMatrices[i],format='csr')
    return matri
#arrays with parameters to build the unitary matrix
def HadamardArraySP(i, k):#i_th wire on totally k wires, start from 0...
    # this should apply Hadamard to wire i out of k wires
    size=2**(k+1)
    #myMatrix=np.zeros((size,size),dtype=complex)
    matrilist = []
    for j in range(k+1):#when j==i tensorproduct with hadamard otherwise by I
        if i != j:
            matrilist.append(scipy.sparse.csr_matrix(scipy.sparse.identity(2,dtype='complex')))
        else:
            matrilist.append(scipy.sparse.csr_matrix(Hadamard()))
    myMatrix = tensorMeSP(matrilist)
    return myMatrix
def PhaseArraySP(i,theta,k,):#i_th wire on totally k wires, start from 0...
    # this should apply Hadamard to wire i out of k wires
    size=2**(k+1)
    matrilist = []
    for j in range(k+1):#when j==i tensorproduct with phase otherwise by I
        if i != j:
            matrilist.append(scipy.sparse.csr_matrix(scipy.sparse.identity(2,dtype='complex')))
        else:
            matrilist.append(scipy.sparse.csr_matrix(Phase(theta)))
    myMatrix = tensorMeSP(matrilist)
    return myMatrix
def CNOTArraySP(ctrl,other,k):#control, other wire on totally k wires, start from 0...
    firstwire = min(ctrl,other)
    matrilist = []
    j = 0
    while j < (k+1):#when j == smaller of the two tensorproduct with cnot, and i += 1, otherwise by I
        if j != firstwire:
            matrilist.append(scipy.sparse.csr_matrix(scipy.sparse.identity(2,dtype='complex')))
            j += 1
        else:
            matrilist.append(scipy.sparse.csr_matrix(CNOT(ctrl,other)))
            j += 2
    myMatrix = tensorMeSP(matrilist)
    return myMatrix
def UnitaryOperationSP(nWires,myInput,vector):#The operation!, and retain 0 or 1 for whether MEASURE
    GateList = []
    matrices = []
    measure = 0
    if myInput[-1] == ['MEASURE']:
        measure = 1
        del myInput[-1]
    for i in range(len(myInput)):
        GateList.append(myInput[i])
    for gate in GateList:
        if gate[0] == 'H':
            matrices.append(HadamardArraySP(int(gate[1]),nWires-1))
        elif gate[0] == 'P':
            matrices.append(PhaseArraySP(int(gate[1]),float(gate[2]),nWires-1))
        elif gate[0] == 'CNOT':
            matrices.append(CNOTArraySP(int(gate[1]),int(gate[2]),nWires-1))
    for matrix in matrices:
        vector = matrix @ vector
    return vector,measure

#3 functions of simulator 1
def Simulator1a(circuitfile):
    UnitaryInfo = UnitaryFromFile(circuitfile)
    Matrix = UnitaryInfo[0][0]
    measure = UnitaryInfo[0][1]
    nwires = UnitaryInfo[2]
    vector = UnitaryInfo[1]
    if len(vector) == 0:
        #print('Using the default basis...')
        vector = np.zeros(2**nwires,dtype=complex)
        vector[0] = 1+0.j
    newvector = np.dot(Matrix,vector)
    newstate = VecToState(newvector)
    if measure == 0:
        return newstate,'Null'
    else:
        Probabilities = MeasureProb(newstate)
        return newstate,Probabilities,ProbToOut(Probabilities)
    pass
def Simulator1b(circuitfile):
    circuit = Read_decompose(circuitfile)
    nwires = circuit[0]
    Input = circuit[1]
    vector = circuit[2]
    if len(vector) == 0:
        #print('Using the default basis...')
        vector = np.zeros(2**nwires,dtype=complex)
        vector[0] = 1+0.j
    if len(vector) == 2**nwires:
        result = UnitaryOperation(nwires,Input,vector)
        newvector = result[0]
        measure = result[1]
        newstate = VecToState(newvector)
        if measure == 0:
            return newstate,'Null'
        else:
            Probabilities = MeasureProb(newstate)
            return newstate,Probabilities,ProbToOut(Probabilities)
    pass
def Simulator1c(circuitfile):
    circuit = Read_decompose(circuitfile)
    nwires = circuit[0]
    Input = circuit[1]
    vector = circuit[2]
    if len(vector) == 0:
        #print('Using the default basis...')
        vector = np.zeros(2**nwires,dtype=complex)
        vector[0] = 1+0.j
        #print(vector)
    if len(vector) == 2**nwires:
        vsparse = scipy.sparse.csr_matrix(vector)
        result = UnitaryOperationSP(nwires,Input,vector)
        newvector = result[0]
        measure = result[1]
        newstate = VecToState(newvector)
        if measure == 0:
            return newstate,'Null'
        else:
            Probabilities = MeasureProb(newstate)
            return newstate,Probabilities,ProbToOut(Probabilities)
    pass

#write circuit and state file to directory for speed test 
def QubitGenerator(wires,filename=''):
    length = 2**wires
    f = open(f'project1/generatedtest/tempstate_{wires}.txt', 'w')
    p = 1
    n = 0
    pterm = []
    while n < length:
        if n == (length-1):
            dice = p
        else:
            dice = np.absolute(np.random.normal((p/(length-n)),(1/length/(5)),1)[0])
            if p < dice:
                dice = 0
            else: pass
        ratio = np.random.uniform()
        real = dice * ratio
        imag = dice - real
        realroot = np.sqrt(real)
        imagroot = np.sqrt(imag)
        if np.random.uniform() > 0.5:
            realroot *= -1
        if np.random.uniform() > 0.5:
            imagroot *= -1
        pterm.append([realroot,imagroot])
        if n == (length-1):
            f.write(str(realroot)+' '+str(imagroot))
        else: f.write(str(realroot)+' '+str(imagroot)+'\n')
        p -= dice
        n += 1
    f.close()
    return pterm
def CircuitGenerator(wires,ngate=100,filename=''):
    f = open(f'project1/generatedtest/tempcircuit_{wires}.circuit', 'w')
    f.write(str(wires)+'\n')
    for i in range(100):
        gatestr = ''
        gate = np.random.randint(1,4)
        wire1 = np.random.randint(0,wires-1)
        if gate == 1: #'H gate'
            gatestr += ('H '+ str(wire1))
        if gate == 2: #'Phase gate'
            gatestr += ('P '+ str(wire1) + ' ' + str(np.random.uniform()*2*np.pi))
        if gate == 3: #'CNOT gate'
            flipped = np.random.randint(0,2)
            wire2 = wire1 + 1
            if flipped == 0:
                gatestr += ('CNOT '+ str(wire1) + ' ' + str(wire2))
            else: gatestr += ('CNOT '+ str(wire2) + ' ' + str(wire1))
        f.write(gatestr + '\n')
    f.write('MEASURE')
    pass

#test the speed of simulators using different circuit and states (same for each choice of nwires)
def TestTime(simulator,nwires):#simulator1 is a string like '1a'
    time0 = time.time()
    if simulator == '1a':
        data = Simulator1a(f"project1/generatedtest/tempcircuit_{nwires}.circuit")
    if simulator == '1b':
        data = Simulator1b(f"project1/generatedtest/tempcircuit_{nwires}.circuit")
    if simulator == '1c':
        data = Simulator1c(f"project1/generatedtest/tempcircuit_{nwires}.circuit")
    if simulator == '2':
        data = Simulator2(f"project1/generatedtest/tempcircuit_{nwires}.circuit")
    time1 = time.time() - time0
    return time1

#begin simulator 2
def RemoveDuplicates(myState):
    newstate = list([0,0] for i in range(2**len(myState[0][1])))
    for element in myState:
        index = int(element[1],2)
        newstate[index][0] += element[0]
        newstate[index][1] = element[1]
    newstate[:] = (value for value in newstate if value != [0,0])
    return newstate
def Hadamard2(inwire,nwires,instate):
    newstate = []
    for element in instate:
        oldstr = element[1]
        oldvalue = element[0]
        change = oldstr[inwire]
        newval = []
        newchange = ['0',(-1)**int(change)]
        for stuff in newchange:
            if stuff == -1:
                newpair = [-1,'1']
            else: newpair = [1,str(stuff)]
            if inwire == 0:
                newpair[1] = newpair[1] + oldstr[1:]
            else:
                newpair[1] = oldstr[:inwire] + newpair[1] + oldstr[inwire+1:]
            newpair[0] *= (oldvalue/np.sqrt(2))
            newval.append(newpair)
        for elem in newval:
            newstate.append(elem)
    return RemoveDuplicates(newstate)
def Phase2(inwire,theta,nwires,instate):
    newstate = []
    for element in instate:
        oldstr = element[1]
        oldvalue = element[0]
        change = oldstr[inwire]
        if change == '0':
            newpair = [1*oldvalue,'0']
        else:
            newpair = [np.exp(theta*1.j)*oldvalue,'1']
        if inwire == 0:
            newpair[1] = newpair[1] + oldstr[1:]
        else:
            newpair[1] = oldstr[:inwire] + newpair[1] + oldstr[inwire+1:]
        newstate.append(newpair)
    return newstate
def CNOT2(ctrlwire,otherwire,nwires,instate):
    newstate = []
    for element in instate:
        oldstr = element[1]
        oldvalue = element[0]
        ctrl = oldstr[ctrlwire]
        newpair = []
        if ctrl == '1': 
            change = oldstr[otherwire]
            if change == '1':
                newindex = '0'
            else: newindex = '1'
            newpair = [1*oldvalue,newindex]
            if otherwire == 0:
                newpair[1] = newpair[1] + oldstr[1:]
            else:
                newpair[1] = oldstr[:otherwire] + newpair[1] + oldstr[otherwire+1:]
            newstate.append(newpair)
        else:
            newstate.append(element)
    return newstate
def Func2(inwire,rangewire,nwires,x,N,instate):
    newstate = []
    for element in instate:
        oldstr = element[1]
        #print(oldstr)
        #stringrange = np.linspace(inwire,inwire+rangewire-1,rangewire)
        if len(oldstr) == rangewire:
            intmstring = oldstr
        else:
            intmstring = oldstr[inwire:inwire+rangewire]
        size = 2**rangewire
        j = BinaryToInt(intmstring)
        if j < N:
            k = (j*x)%N
        else:
            k = j
        newpart = str(IntToBinary(k))
        while len(newpart) < len(intmstring):
            newpart = '0' + newpart
        if inwire == 0:
            newstr = newpart + oldstr[rangewire:]
        else:
            newstr = oldstr[:inwire] + newpart + oldstr[inwire+rangewire:]
        #print(size,intmstring,j,k,newpart,newstr)
        newstate.append([element[0],newstr])
    return newstate
def CFunc2(ctrlwire,inwire,rangewire,nwires,x,N,instate):
    newstate = []
    for element in instate:
        oldstr = element[1]
        ctrl = oldstr[ctrlwire]
        if ctrl == '1':
            intmstring = oldstr[inwire:inwire+rangewire]
            size = 2**rangewire
            j = BinaryToInt(intmstring)
            if j < N:
                k = (j*x)%N
            else:
                k = j
            newpart = str(IntToBinary(k))
            while len(newpart) < len(intmstring):
                newpart = '0' + newpart
            if inwire == 0:
                newstr = newpart + oldstr[rangewire:]
            else:
                newstr = oldstr[:inwire] + newpart + oldstr[inwire+rangewire:]
            #print(size,intmstring,j,k,newpart,newstr)
        else:
            newstr = oldstr
        newstate.append([element[0],newstr])
    return newstate
def Simulator2(circuitfile,debugging=0):
    #global T0,T1,T2,T3
    circuit = Read_decompose(circuitfile)
    nwires = circuit[0]
    Input = circuit[1]
    vector = circuit[2]
    if len(vector) == 0:
        #print('Using the default basis...')
        vector = np.zeros(2**nwires,dtype=complex)
        vector[0] = 1+0.j
        
    if len(vector) == 2**nwires:
        resultstate = VecToState(vector)
        if debugging > 1:
            print(resultstate)
        measure = 0
        if Input[-1][0] == 'MEASURE':
            measure = 1
        for elem in Input:
            if elem[0] == 'H':
                #tx = time.time()
                resultstate = Hadamard2(int(elem[1]),nwires,resultstate)
                #T0 += time.time()-tx
            elif elem[0] == 'P':
                #tx = time.time()
                resultstate = Phase2(int(elem[1]),eval(str(elem[2])),nwires,resultstate)
                #T1 += time.time()-tx
            elif elem[0] == 'CNOT':
                #tx = time.time()
                resultstate = CNOT2(int(elem[1]),int(elem[2]),nwires,resultstate)
                #T2 += time.time()-tx
            elif elem[0] == 'FUNC': #the classical function!!!!!
                resultstate = Func2(int(elem[1]),int(elem[2]),nwires,int(elem[4]),int(elem[5]),resultstate)
            elif elem[0] == 'CFUNC': #the classical function!!!!!
                #tx = time.time()
                resultstate = CFunc2(int(elem[1]),int(elem[2]),int(elem[3]),nwires,int(elem[5]),int(elem[6]),resultstate)
                #T3 += time.time()-tx
            else:pass
        resultstate.sort(key=lambda x:x[1]) #reference: https://www.pythonpool.com/python-sort-list-of-lists/
        resultstate = VecToState(StateToVec(resultstate))
        if measure == 0:
            return resultstate,'Null'
        else:
            Probabilities = MeasureProb(resultstate)
            Out = ProbToOut(Probabilities)
            if debugging >= 1:
                return resultstate,Probabilities,Out
            else:
                return Out
    pass

#simulate measurements from simulators on a histogram
def Histogram(prob,n):
    num = []
    for char in np.linspace(0,(len(prob)-1),len(prob)):
        num.append(str(IntToBinary(int(char))))
    count = np.zeros(len(prob))
    for i in range(n):
        count[BinaryToInt(ProbToOut(prob))] += 1
    plt.figure(figsize=(7,7))
    plt.bar(num,count,width=1/32*len(count),align='center')
    plt.xticks(rotation = 90)
    plt.xlabel('state')
    plt.ylabel(f'count,total={n}')
    plt.show()
    pass

#phase estimation
def PEst(filename,BASIS,N,topwire,bottomwire=0,size=5): #no .circuit here in filename
    #https://pythonexamples.org/python-replace-string-in-file/#:~:text=To%20replace%20a%20string%20in%20File%20using%20Python,,file.%204%20Close%20both%20input%20and%20output%20files.
    angles = np.linspace(0,2*np.pi,N)
    result = []
    resultb = []
    for i in range(N):
        angle = angles[i]
        #input file
        fin = open(f'{filename}.circuit', "rt")
        fout = open(f'{filename}_tp.circuit', "wt")
        for line in fin:
            line = line.replace('X', str(angle))
            line = line.replace('Y', str(BASIS))
            fout.write(line)
        fin.close()
        fout.close()
        output = BinaryToInt(Simulator2(f'{filename}_tp.circuit',0)[:-1*(bottomwire)])
        result.append(output/(2**(topwire)))
        p = Simulator2(f'{filename}_tp.circuit',1)[1]
        maxp = 0
        maxs = ''
        for elem in p:
            if elem[1] > maxp:
                maxs = elem[0]
                maxp = elem[1]
        resultb.append(BinaryToInt(maxs[:-1*(bottomwire)])/(2**(topwire)))
    plt.figure(figsize=(size,5))
    #plt.ylim((0,1))
    plt.xlabel('phi/(2*pi)')
    plt.ylabel('theta j')
    plt.plot(angles/(2*np.pi),resultb)
    plt.savefig(f'project1/files/{topwire}_{BASIS[1:-1]}_count={N} phases.png')
    pass
def PEstHistogram(filename,BASIS,N,topwire,bottomwire=1,theta=(2*np.pi*0.1432394487827058)):
    result = []
    for i in range(N):
        fin = open(f'{filename}.circuit', "rt")
        fout = open(f'{filename}_tp.circuit', "wt")
        for line in fin:
            line = line.replace('X', str(theta))
            line = line.replace('Y', str(BASIS))
            fout.write(line)
        fin.close()
        fout.close()
        R = Simulator2(f'{filename}_tp.circuit',0)[:-1*(bottomwire)]
        #print(R)
        result.append(BinaryToInt(R)/(2**topwire))
    x = np.linspace(0,(1-(0.5**topwire))*2*np.pi,2**topwire)
    
    xp = x/2/np.pi
    #print(xp)
    y = np.zeros(2**topwire)
    for i in range(len(result)):
        k = result[i]
        y[int((2**topwire)*k)] += 1
    #print(y)
    y/=N
    maxy = 0
    for m in y:
        maxy = max(maxy,m)
    xt = [theta/2/np.pi,theta/2/np.pi]
    yt = [0,maxy]
    plt.figure()
    plt.xlim((-0.05,1.05))
    plt.bar(xp,y,width=0.03125/2,align='center')
    plt.plot(xt,yt,'k-')
    plt.xlabel(f'output,theta/2pi={theta/2/np.pi}')
    plt.ylabel('probability')
    #plt.show()
    plt.savefig(f'project1/files/{topwire}_{BASIS[1:-1]}_{theta}_count={N} histogram.png')
    pass
def INVQFT(nwires):#inversed QFT of n wires with SWAPs after. Add this after the input run by simulator
    #after swapping and then build is easier for the indexing!!!
    #so, QFT template is the reverse of wikipedia one
    
    gatelist = []
    wire = 0 #this number is index (start from 0), not the count of wires.
    #make the inversed QFT gates running backward (upside-down >leftright< of wikipedia thing, similar to 2 wire example)
    while wire < nwires:
        remain = 0
        for i in range(wire):
            if remain <= wire:
                gatelist.append(['CPHASE',wire-remain-1,wire,((0.5)**(wire-remain))*np.pi*-1])
                #*-1--->negative pi/2**X is the inverse of CP matrix
                remain += 1
        gatelist.append(['H',wire])
        wire += 1
    Twire = 0
    Bwire = nwires-1
    while Twire < Bwire:
        if Twire != Bwire:
            gatelist.append(['SWAP',Twire,Bwire])
        Twire += 1
        Bwire -= 1

    #NOW SWAP all wires at the end:
    
    return gatelist
def INVQFT2(nwires):#swap before and QFT after, seems like more accurate
    
    gatelist = []

    Twire = 0
    Bwire = nwires-1
    while Twire < Bwire:
        if Twire != Bwire:
            gatelist.append(['SWAP',Twire,Bwire])
        Twire += 1
        Bwire -= 1

    wire = nwires-1
    while wire >= 0:
        remain = wire
        gatelist.append(['H',wire])
        for i in range(wire):
            while remain < nwires:

                gatelist.append(['CPHASE',wire-1,nwires-remain+wire-1,((0.5)**(nwires-remain))*np.pi*-1])
                #*-1--->negative pi/2**X is the inverse of CP matrix
                remain += 1
        
        wire -= 1
    #NOW SWAP all wires at the end:
    return gatelist

def AddInvQFT(filename,nwires):
    fin = open(f'{filename}.circuit', "rt")
    fout = open(f'{filename}_qft.circuit', "wt")
    measure = 0
    for line in fin:
        if 'MEASURE' not in line:
	        fout.write(line)
        else:
            measure = 1
    QFT = INVQFT2(nwires)
    for gate in QFT:
        line = ''
        for element in gate:
            line += (str(element)+' ')
        line += '\n'
        fout.write(line)
    if measure == 1:
        fout.write('MEASURE')
    fin.close()
    fout.close()
    pass

#construct quantum Shor's algorithm
def newcircuit(x,N,debug=0,f=''):
    if f == 'save':
        f = open(f'project1/shors_{x}_{N}.circuit', 'w')
    else:
        f = open(f'project1/shorstemp.circuit', 'w')
    n = int(np.ceil((np.log10(N)/np.log10(2))))
    #print(n)
    wires = (2*n + 1)
    f.write(str(wires)+'\n')
    if debug == 0:
        basis = '1'
        while len(basis) < wires:
            basis = '0' + basis
        basis = ('|' + basis + '>')
        f.write(f'INITSTATE BASIS {basis}\n')
    else:
        f.write(f'INITSTATE BASIS Y\n')
    topwire = n+1
    bottomwire = n
    #generate Hadamards
    for i in range(topwire):
        f.write(f'H {i}\n')
    for j in range(topwire):
        rep = 2**j
        f.write(f'CFUNC {topwire-1-j} {topwire} {n} xyModN {x**(2**(j))} {N}\n')
        #for r in range(rep):
            #f.write(f'CFUNC {topwire-1-j} {topwire} {n} xyModN {x} {N}\n')
    f.write('MEASURE')
    f.close()
    if f == 'save':
        AddInvQFT(f'project1/shors_{x}_{N}',topwire)
    else:
        AddInvQFT(f'project1/shorstemp',topwire)
    return wires,topwire,bottomwire
def Shors(x,N):
    wireconfig = newcircuit(x,N)
    #d2 = Simulator2('project1/shorstemp_qft.circuit',1)
    #print(d2[0],x)
    R = Simulator2('project1/shorstemp_qft.circuit')[:-1*(wireconfig[2])]
    decimal = BinaryToInt(R)/(2**(wireconfig[1]))
    return decimal
def FACTOR(N):
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
                value = Shors(x,N)
                r = Fraction(value).limit_denominator(100).denominator
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

#predefined states
myState=[(np.sqrt(0.1), '00'),(np.sqrt(0.4), '01') , (-np.sqrt(0.5), '11' )]
myState2=[ (np.sqrt(0.1)*1.j, '101'), (np.sqrt(0.5), '000') , (-np.sqrt(0.4), '010' ) ]

print('---Start Run---')
#put commands here
print(FACTOR(21))

print('---End Run---')
