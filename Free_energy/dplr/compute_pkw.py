import numpy as np
import matplotlib.pyplot as plt

c0=6.023E23/1.E27 #Standard concentration in atoms/A3

def load_fes(sys):
  f=np.loadtxt("{}/fes.dat".format(sys)).T
  return f

def integrate(F,S,Sc):
  intF=4.E0*np.pi*np.trapz(S[S<Sc]**2.E0*np.exp(-F[S<Sc]),S[S<Sc])*(S[-1]-S[-2])
  return intF

def compute_pkw(f,Ssep):
  pHinv=c0*integrate(f[1],f[0],Ssep)
  pH=np.log10(pHinv)
  return 2*pH
  
def main():
  sys=[64,128,256,512,1024]
  pKw=np.zeros(len(sys))
  l=4.0
  for s in range(len(sys)):
    fes=load_fes(sys[s])
    fes[1]-=np.mean(fes[1][-5:])
    pKw[s]=compute_pkw(fes,l) 

def plot_curves():
  plt.xlabel('Number of H$_2$O molecules',fontsize=9)
  plt.ylabel('pKw',fontsize=9)
  plt.legend()
  plt.tight_layout()
  plt.show()

fig,ax=plt.subplots(1,1,figsize=(3.3,3.3))
main()
plot_curves()
