![image](LGCMC-logo.png)

# LGCMC Code Repository

Welcome to **LGCMC** - A FORTRAN package for Lattice Grand Canonical Monte Carlo Simulation for *ab initio* prediction of pure and mixed gas isotherms.

## Introduction

LGCMC employs the Lattice-Gas Hamiltonian:

```math
 H(n_{1},n_{2},...,n_{M}) = \sum_{s=1}^S\sum_{c=A}^K\Delta G_{s}^{c}n_{i}^{c} + 
 \frac{1}{2} \sum_{s \neq s'} \sum_{c,c'=A}^K E_{s,s'}^{c,c'} 
n_{s}^{c}n_{s'}^{c'} + ... 
```
to define the free energies of a multi-component system on a surface. 
Here, $s,s'=1,2,...,M$ denote the indices of adsorption sites within the 
simulation box, where $M$ is the maximum number of available sites. 
The indices $c,c'=A,B,...,K$ represent the gas species of a multi-component 
mixture. The Gibbs free energy of adsorption $(G_{s}^{c})$ of component $c$ 
on an isolated site $s$, and lateral interaction energies ($E_{s,s'}^{c,c'}$) 
of components $c$ and $c'$ adsorbed at sites $s$ and $s'$, respectively, 
can be computed applying accurate quantum chemical methods. The occupancy number
of the $c$-th species at the $s$-th site, $n_{s}^{c}$, is related to the 
total occupancy, $n_{s}$, according to:

```math
\sum_{c=A}^{K} n_{s}^{c} = n_{s}.
```

An adsorption site can be occupied by atmost one gas molecule. Therefore, both 
$n_{s}^{c}$ and $n_{s}$ can have only two allowed values: 0 and 1. 

The number of gas molecules of species $c$ at a specific configuration $i$ is:
```math
N_{c,i}=\sum_{s=1}^{M} n_{s}^{c}
``` 
LGCMC generates samples configurations: $i=1,2,3,...,I$, using a grand canonical 
Monte Carlo simulation where chemical potential of each species is same in both gas and solid phase, *i.e.*, $\mu_{c}^{ads} = \mu_{c}^{gas}$. 
The condition of materials equilibrium dictates that the gas phase chemical 
potential, $\mu_{c}^{gas} = f(y,P,T)$, which is a thermodynamic function 
of pressure $(P)$, temperature $(T)$, and gas phase composition $(y)$, 
uniquely defines the adsorbed phase chemical potential as well. 
The chemical potential can be computed employing a suitable equation of state or
from the experimentally measured fugacity coefficients: 
$\phi_{c}=e^{\frac{\mu_{c}}{RT}}$ with $R$ denoting the universal gas constant.  The following equation of states: (i) ideal, (ii) Van der Waals, (iii) Redlich-Kwong, (iv) Soave-Redlich-Kwong, (v) Peng-Robinson (with 1978 and 1980's 
modifications), (vi) Peng-Robinson-Gasem equation of states are implemented. 
It is also possible to supply the fugacity coefficients (if known from
 experiments) to LGCMC.

After generating $I$ configurations, the surface coverage of component $c$ is 
defined as:

```math
\theta_{c}=\frac{\langle N_{c} \rangle}{M} = \frac{1}{MI}\sum_{i=1,I}N_{c,i} 
``` 

 
## Installation
**Requirements:** Intel fortran compiler.
To compile,
(1) Change install.sh to an executable by chmod 755 install.sh.
(2) Load intel fortran compiler by "module load <compiler version>".
(3) Run install.sh using ./install.sh

Manual is comming soon.
