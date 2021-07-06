import numpy as np
import matplotlib.pyplot as plt

G = 4*np.pi**2 # Gravitational constant
c = 63197.8 # Speed of light in units of au/yr

# Masses  
m1 = 10. # Solar masses
m2 = 7.  # Solar masses

# Radii
r1 = 2*G*m1/c**2 # Schwarzschild radius for m1
r2 = 2*G*m2/c**2 # Schwarzschild radius for m2
r_crit = r1 + r2 # Separation distance for fusion/collision

M = m1 + m2 # Total mass (solar masses)
GM = G*M # gravitational parameter


# Initial Values
E = -70. # energy
L = 50. # angular momentum
omega = 0. #np.pi/3 # argument of the pericenter



# Orbital parameters
def semilatus(L):
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

def eccentricity(E, L):
  '''
  Returns the eccentricity of the orbit from the 
  energy E and the angular momentum L
  '''
  return np.sqrt(1 + 2*E*L**2/(GM**2))


def phi(L, r, ang0):
  '''
  Returns the true anomaly obtained from E, L and r
  '''
  ang1 = ang0 + dt*L/r**2
  return ang1

def ellipse(p, ecc, f):
  '''
  Returns the coordinates of the point and 
  the radius at angle f
  '''
  r = p/(1-ecc*np.cos(f - omega))
  x = r*np.cos(f)
  y = r*np.sin(f)
  return x, y, r


# Grid definition
n = 10000 # number of nodes
time = np.linspace(0., 10., n)
dt = time[1] - time[0] # timestep
position = np.zeros([n,3])
angle = np.zeros(n)

# Initial condition in the grid
angle[0] = 0.
position[0] = ellipse(semilatus(L), eccentricity(E,L), angle[0])



# Main Loop
for i in range(n-1):
  ecc = eccentricity(E,L)
  p = semilatus(L)
  # Check for ellipse's eccentricity
  if ecc<0 or ecc>=1:
    print('The eccentricity is not that of an ellipse.')
    break
  angle[i+1] = phi(L, position[i,2], angle[i])
  position[i+1] = ellipse(p, ecc, angle[i+1])
  # First estimate for fusion of the binary system
  if position[i+1,2]<1E5*r_crit:
    print('The binary system has merged at time {:.0f}'.format(i))
    break
  #Energy and Angular Momentum lose model 
  L = L-0.005
  E = E-0.005


#print(r_crit)

#print(ecc(E,L), a(E))
#print(np.sqrt(2*E+2*GM))
#print(dt)
#print(E,L)

u= 8270
w = 8300
plt.figure(figsize=(10,7))
plt.plot(position[:,0],position[:,1])
#plt.plot(position[u:w,0],position[u:w,1])
plt.axhline(color='black',alpha=0.3)
plt.axvline(color='black',alpha=0.3)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.show()