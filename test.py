import numpy as np
string = '|000>'
print(string[1:-1])
list1 = []
print(len(list1))
#INITSTATE BASIS |100>

#stt= input(f'statefile is: ')
#print(stt)
print((eval(str(np.pi))))

string1 = 'algbnelwgjsd'
string1 = string1.replace('a','X')
string1 = string1.replace('b','Y')
print(string1)

N = 4
x = np.linspace(0,1-(0.5**N),2**N)
print(x)
print(2*np.pi*0.1432394487827058)
print(int('1111111',2))
print(np.sqrt(0.3),np.sqrt(0.7))

print(np.kron(((1,0),(0,4)),np.identity(2)))
print(np.kron(np.identity(2),((1,0),(0,4))))