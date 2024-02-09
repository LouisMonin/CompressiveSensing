# CompressiveSensing
DCT function to generate the matrix of a transformed into discrete cosine of size N

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
def DCT(N):
	C=np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			if j==0:
				C[i,j]=1/np.sqrt(N)
			else:
				C[i,j]=np.sqrt(2/N)*np.cos(j*(2*i+1)*np.pi/(2*N))
	return C

#Pour vérifier que la matrice donnée par DCT(N) est bien orthogonale, on la multiplie par sa transposée et on vérifie que l'on obtient bien l'identité.

C=DCT(10)
print(np.dot(C,C.T))

print("")
print("")

t=np.arange(0,499.5/430,1/430)
print(len(t))

f0=50
fp=100
s=np.sin(2*np.pi*f0*t)*np.cos(2*np.pi*fp*t)
print(len(s))
plt.plot(s)
plt.show()

print("")
print("")

C=DCT(500)
alpha=(np.dot(C.T,s))
plt.plot(alpha)
plt.show()

print("")
print("")

def DFT(N):
	C=np.zeros((N,N),dtype=complex)
	for p in range(N):
		for k in range(N):
			C[p,k]=np.sqrt(1/N)*np.exp(1j*2*p*k*np.pi/(N))
	return C

print(DFT(3))
print("")
print("")

F=DFT(3)
I=np.dot(F,np.conj(F.T))
print(I)

print("")
print("")

F=DFT(500)
alpha2=np.dot(np.conj(F.T),s)
plt.plot(alpha2)
plt.show()
plt.plot(np.abs(alpha2))
plt.show()
