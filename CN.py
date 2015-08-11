"""
	This program solves the heat equation
		u_t = u_xx
	with dirichlet boundary condition
		u(0,t) = u(1,t) = 0
	with the Initial Conditions
		u(x,0) = 5*cos( pi*x )**2
	over the domain x = [0, 1]
 
	The program solves the heat equation using a finite difference method
	where we use a center difference method in space and Crank-Nicolson 
	in time.
"""
 
import scipy as sc
import scipy.sparse as sparse
import scipy.sparse.linalg
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
 
import mpl_toolkits.mplot3d


# Number of internal points
N = 200
 
# Calculate Spatial Step-Size
h = 1/(N+1.0)
 
# Create Temporal Step-Size, TFinal, Number of Time-Steps
k = h/2
TFinal = 0.3
NumOfTimeSteps = int(TFinal/k)
print "num of time",NumOfTimeSteps 
# Create grid-points on x axis
x = np.linspace(0,1,N+2)
x = x[1:-1]
 
# Initial Conditions
u = np.transpose(np.mat(5*np.cos(np.pi*x)**2))
 
# Second-Derivative Matrix
data = np.ones((3, N))
data[1] = -2*data[1]
#boundary
data[2,1]=2*data[2,1]
data[0,N-2]=2*data[0,N-2]
diags = [-1,0,1]
D2 = sparse.spdiags(data,diags,N,N)/(h**2)
 
# Identity Matrix
I = sparse.identity(N)
 
# Data for each time-step
data = []
 
for i in range(NumOfTimeSteps):
	# Solve the System: (I - k/2*D2) u_new = (I + k/2*D2)*u_old
	A = (I -k/2*D2)
	b = ( I + k/2*D2 )*u
	u = np.transpose(np.mat( sparse.linalg.spsolve( A,  b ) ))
 
	# Save Data
	data.append(u)
 

print "u=",u
print "shape of data",np.shape(data) 
aa=np.shape(data)
data=np.reshape(data,(aa[0],aa[1]))
print "shape of data",np.shape(data) 
# Define the Frame Speed and Movie Length
plt.plot(x, data[1],'ro' )
plt.plot(x, data[10],'bo' )
plt.plot(x, data[30],'y^' )
plt.plot(x, data[40],'k^' )
#plt.imshow(data)

plt.show()
 
 
