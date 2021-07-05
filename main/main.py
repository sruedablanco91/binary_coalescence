import numpy as np
import matplotlib.pyplot as plt

class ellipse:
  def __init__(self, a=1., e=0.0167, omega=0.):
    if e >=1. or e<0.:
      self.is_ellipse = False
      self.mes = f'The eccentricity {self.e:.5f} does not correspond to an ellipse!'
      print(self.mes)
    else:
      self.is_ellipse = True
      self.a = a
      self.e = e
      self.omega = omega
      self.semilatus = self.a*(1-self.e**2)
      self.r_min = self.a*(1-self.e)
      self.r_max = self.a*(1+self.e)
  
  def plot(self, color='crimson'):
    if self.is_ellipse == True:
      boundary = 2.*self.a
      f = np.linspace(0, 2*np.pi, 1000)
      x = self.semilatus*np.cos(f)/(1+self.e*np.cos(f-self.omega))
      y = self.a*(1-self.e**2)*np.sin(f)/(1+self.e*np.cos(f-self.omega))

      xsemaxis = np.linspace(-self.r_max*np.cos(self.omega), 
                             self.r_min*np.cos(self.omega), 100)
      ysemiaxis = np.tan(self.omega)*xsemaxis
      
      plt.figure(figsize=(7,7))
      plt.plot(x,y, color=color)
      plt.plot(xsemaxis,ysemiaxis, '--', color='grey')
      plt.axhline(color='black',alpha=0.3)
      plt.axvline(color='black',alpha=0.3)
      plt.xlim(-boundary, boundary)
      plt.ylim(-boundary, boundary)
      plt.xlabel(f'$x$')
      plt.ylabel(f'$y$')
      plt.show()
    else:
      return print(self.mes)
  
