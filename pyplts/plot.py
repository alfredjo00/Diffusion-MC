import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

sns.set('talk')

## Load Data

data1 = np.loadtxt("data/task1_E-T_N.txt")

data2 = np.loadtxt("data/task2_E-T_N.txt")

data3a = np.loadtxt("data/task3a_E-T_N.txt")

data3b = np.loadtxt("data/task3b_E-T_N.txt")

def normalize(x, y):
    norm = np.trapz(y, x=x)
    return y/norm

################################################################

def plot_walkers():
    N1 = data1[:,1]
    N2 = data2[:,1]
    N3a = data3a[:,1]
    N3b = data3b[:,1]
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(N1, label=r'Walkers $N$', linewidth=1)
    ax.legend()
    
    ax.set_ylabel(r"$N$ Walkers")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.02$)")
    fig.tight_layout()
    fig.savefig("graphs/task1_N_walkers.jpg", dpi=120)   
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(N2, label=r'Walkers $N$', linewidth=1)
    ax.legend()
    
    ax.set_ylabel(r"$N$ Walkers")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.01$)")
    ax.set_title("DMC")
    fig.tight_layout()
    fig.savefig("graphs/task2_N_walkers.jpg", dpi=120)   
    
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(N3a, label=r'Walkers $N$', linewidth=1)
    ax.plot([0, len(N3a)], [1000,1000], '--k', label='$N=1000$')
    ax.legend()
    
    ax.set_ylabel(r"$N$ Walkers")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.1$)")
    ax.set_title("DMC with importance sampling")
    fig.tight_layout()
    fig.savefig("graphs/task3a_N_walkers.jpg", dpi=120)   
    
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(N3b, label=r'Walkers $N$', linewidth=1)
    ax.plot([0, len(N3b)], [1000,1000], '--k', label='$N=1000$')
    ax.legend()
    
    ax.set_ylabel(r"$N$ Walkers")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.1$)")
    ax.set_title("DMC with importance sampling, 2nd order drift")
    fig.tight_layout()
    fig.savefig("graphs/task3b_N_walkers.jpg", dpi=120)   
    
    
################################################################

def plot_walkers_dist():
    dist = np.loadtxt("data/task1_walkers.txt")
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.hist(dist, bins=30, label=r'Walkers $N$')
    ax.legend()
    
    fig.tight_layout()
    fig.savefig("graphs/task1_walkers_dist.jpg", dpi=120)   
    
    
################################################################

def plot_E():
    n_skip = 2000
    E1 = data1[:,0]
    mean = np.mean(E1[n_skip:])
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(E1, label=r'$E_T$', linewidth=1)
    ax.plot([0, len(E1)], [mean, mean], '--k', label=fr'$\mu(E_T)={mean:.12f} \, E_H$')
    
    ax.set_ylabel(r"Energy $E_T$ [$E_H$]")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.02$)")
    
    ax.legend()
    
    fig.tight_layout()
    fig.savefig("graphs/task1_E_T.jpg", dpi=120)   
    
    
    E2 = data2[:,0]
    mean = np.mean(E2[n_skip:])
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    
    ax.set_ylabel(r"Energy $E_T$ [$E_H$]")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.01$)")
    
    ax.plot(E2, label=r'$E_T$', linewidth=1)
    ax.plot([0, len(E2)], [mean, mean], '--k', label=fr'$\mu(E_T)={mean:.12f} \, E_H$')
    ax.legend()
    
    ax.set_title("DMC without importance sampling")
    fig.tight_layout()
    fig.savefig("graphs/task2_E_T.jpg", dpi=120)   
        
    
    E3a = data3a[:,0]
    mean = np.mean(E3a[n_skip:])
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    
    ax.set_ylabel(r"Energy $E_T$ [$E_H$]")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.1$)")
    
    ax.plot(E3a, label=r'$E_T$', linewidth=1)
    ax.plot([0, len(E3a)], [mean, mean], '--k', label=fr'$\mu(E_T)={mean:.12f} \, E_H$')
    ax.legend()
    
    ax.set_title("DMC with importance sampling")
    fig.tight_layout()
    fig.savefig("graphs/task3a_E_T.jpg", dpi=120)  
        
    
    E3b = data3b[:,0]
    mean = np.mean(E3b[n_skip:])
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    
    ax.set_ylabel(r"Energy $E_T$ [$E_H$]")
    ax.set_xlabel(r"Time steps ($\Delta \tau = 0.1$)")
    
    ax.plot(E3b, label=r'$E_T$', linewidth=1)
    ax.plot([0, len(E3b)], [mean, mean], '--k', label=fr'$\mu(E_T)={mean:.12f} \, E_H$')
    ax.legend()
    
    ax.set_title("DMC with importance sampling, 2nd order drift")
    fig.tight_layout()
    fig.savefig("graphs/task3b_E_T.jpg", dpi=120)   
    
################################################################

def extrapolate_E():
    
    data_e = np.loadtxt("data/task4_energy.txt")
    
    X = np.linspace(0.01, 0.4, len(data_e[:,0]))
        
        
    z1 = np.polyfit(X, data_e[:,0], deg=1)
    z2 = np.polyfit(X, data_e[:,1], deg=2)
    
    p1 = np.poly1d(z1)
    p2 = np.poly1d(z2)
    
    print(p1(0), p2(0))
    
    f1 = [p1(t) for t in np.linspace(0,0.4,50)] 
    
    f2 = [p2(t) for t in np.linspace(0,0.4,50)] 
    
    
    
    fig,ax = plt.subplots(figsize=(8,6))
    
    ax.plot(X, data_e[:,0], '--r', label="1st composition")    
    
    ax.plot(X, data_e[:,1], '-.b', label="2nd composition")
    
    ax.plot(np.linspace(0,0.4,50), f1, '-k', label='Fit of 1st')
    
    ax.plot(np.linspace(0,0.4,50), f2, ':k', label='Fit of 2nd')    
    
    ax.legend()
    
    fig.tight_layout()
    fig.savefig("graphs/energy_extrapolation.jpg", dpi=100)
    
################################################################

if __name__ == '__main__':
    plot_walkers()
    plot_E()
    plot_walkers_dist()
    extrapolate_E()
    