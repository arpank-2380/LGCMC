![image](LGCMC-logo.png)

# LGCMC Code Repository

Welcome to **LGCMC** - A FORTRAN package for Lattice Grand Canonical Monte Carlo Simulation for *ab initio* prediction of pure and mixed gas isotherms.

## Introduction

LGCMC employs Lattice-Gas Hamiltonian:

```math
 H(n_{1},n_{2},...,n_{M}) = \sum_{s=1}^M\sum_{c=A}^K\Delta G_{s}^{c}n_{i}^{c} + 
                           \frac{1}{2}\sum_{s!=s'}\sum_{c,c'=A}^KE_{s,s'}^{c,c'} 
```

Requirements: Intel fortran compiler.
To compile,
(1) Change install.sh to an executable by chmod 755 install.sh.
(2) Load intel fortran compiler by "module load <compiler version>".
(3) Run install.sh using ./install.sh

Manual is comming soon.
