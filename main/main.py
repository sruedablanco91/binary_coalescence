import numpy as np
import matplotlib.pyplot as plt

# Gravitational constant
G = 4*np.pi**2

# Masses  
m1 = 10. # Solar masses
m2 = 7.  # Solar masses

M = m1 + m2 # Total mass (solar masses)
GM = G*M # gravitational parameter


# Initial Values
E = -70. # energy
L = 50. # angular momentum
omega = np.pi/3
x0 = 1. # au
y0 = 0. # au


# Orbital parameters
def p(L):
  '''
  Returns the semilatus rectum of the orbit 
  from the angular momentum L
  '''
  return L**2/GM

def a(E):
  '''
  Returns the semimajor axis of the orbit 
  in au from the energy E
  '''
  return -GM/(2*E)

def ecc(E, L):
  '''
  Returns the eccentricity of the orbit from the 
  energy E and the angular momentum L
  '''
  return np.sqrt(1 + 2*E*L**2/(GM**2))


def phi(E, L, pos0, ang0):
  '''
  Returns the true anomaly obtained from E, L and r
  '''
  r = np.linalg.norm(pos0) 
  ang1 = ang0 + dt*L/r**2
  return ang1

def ellipse(E, L, f):
  '''
  Returns the coordinates of the point at angle f
  '''
  x = p(L)*np.cos(f)/(1-ecc(E,L)*np.cos(f - omega))
  y = p(L)*np.sin(f)/(1-ecc(E,L)*np.cos(f - omega))
  return x, y


# Grid definition
n = 10000 # number of nodes
time = np.linspace(0., 10., n)
dt = time[1] - time[0] # timestep
position = np.zeros([n,2])
angle = np.zeros(n)

# Initial condition in the grid
angle[0] = 0.
position[0] = ellipse(E, L, angle[0])



# Main Loop
for i in range(n-1):
  #L = L-0.005
  #E = E-0.005
  if ecc(E,L)<0 or ecc(E,L)>=1:
    print('The eccentricity is wrong')
    break
  angle[i+1] = phi(E, L, position[i], angle[i])
  position[i+1] = ellipse(E, L, angle[i+1])


position[0:3]

#print(ecc(E,L), a(E))
#print(np.sqrt(2*E+2*GM))
#print(dt)


plt.figure(figsize=(7,7))
plt.plot(position[:,0],position[:,1])
plt.axhline(color='black',alpha=0.3)
plt.axvline(color='black',alpha=0.3)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.show()