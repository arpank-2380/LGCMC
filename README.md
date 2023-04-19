![image](LGCMC-logo.png)

# LGCMC Code Repository

Welcome to **LGCMC** - A FORTRAN package for Lattice Grand Canonical Monte Carlo Simulation for *ab initio* prediction of pure and mixed gas isotherms.

## Introduction

LGCMC employs the Lattice-Gas Hamiltonian:

```math
 H(n_{1},n_{2},...,n_{S}) = \sum_{s=1}^S\sum_{c=A}^K\Delta G_{s}^{c}n_{i}^{c} + 
 \frac{1}{2} \sum_{s \neq s'} \sum_{c,c'=A}^K E_{s,s'}^{c,c'} 
n_{s}^{c}n_{s'}^{c'} + ... 
```
to define the free energies of a multi-component system on a surface. 
Here, $s,s'=1,2,...,S$ denote the indices of adsorption sites within the 
simulation box, where $S$ is the maximum number of available sites. 
The indices $c,c'=A,B,...,K$ represent the gas species of a multi-component 
mixture. The Gibbs free energy of adsorption of component $c$ on an isolated 
site $s$ $(G_{s}^{c})$, and lateral interaction energies ($E_{s,s'}^{c,c'}$) 
of components $c$ and $c'$ adsorbed at sites $s$ and $s'$, respectively, 
can be computed applying accurate quantum chemical methods. The occupancy number
of the $c$-th species at the $s$-th site, $n_{s}^{c}$, is related to the 
total occupancy, $n_{s}$, according to:

```math
\sum_{c=A}^{K} n_{s}^{c} = n_{s}.
```

An adsorption site can be occupied by atmost one gas molecule. Therefore, both 
$n_{s}^{c}$ and $n_{s}$ can have only two allowed values: 0 and 1.

 

Requirements: Intel fortran compiler.
To compile,
(1) Change install.sh to an executable by chmod 755 install.sh.
(2) Load intel fortran compiler by "module load <compiler version>".
(3) Run install.sh using ./install.sh

Manual is comming soon.
