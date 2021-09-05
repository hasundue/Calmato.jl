# Calmato

[![CI](https://github.com/hasundue/Calmato.jl/workflows/CI/badge.svg)](https://github.com/hasundue/Calmato.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/hasundue/Calmato.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hasundue/Calmato.jl)

**Calmato** is an open-source CALPHAD (CALculation of PHAse Diagrams) package in Julia language, which utilizes [**EAGO**](https://github.com/PSORLab/EAGO.jl) as its optimization engine. The package is still in the alpha stage of development and has limited functionality and reliability.

## Example

Read a tdb file:

```julia
using Calmato

db = read_tdb("cuzn_liang.tdb")
```

```terminal
Database Summary:
        Elements: 4
        Functions: 10
        Phases: 8
```

Print contents of the database:

```julia
print(db)
```

```terminal
Database:
        Elements: 4
                /-
                VA
                Cu
                Zn
        Functions: 10
                RTLNP
                GHSERCU
                GHSERZN
                GLIQCU
                GBCCCU
                GHCPCU
                GLIQZN
                GBCCZN
                GFCCZN
                A1BCZ
        Phases: 8
                1: Liquid; (Cu,Zn)1
                        G(Cu;0)
                        G(Zn;0)
                        L(Cu,Zn;0)
                        L(Cu,Zn;1)
                        L(Cu,Zn;2)
                2: Bcc; (Cu,Zn)1
                        G(Cu;0)
                        G(Zn;0)
                        L(Cu,Zn;0)
                        L(Cu,Zn;1)
                        L(Cu,Zn;2)
                3: Bcc_B2; (Cu,Zn)0.5(Cu,Zn)0.5
                        ...
```

Initialize a system (automatically removing electron gas, "/-", and vacancy, "VA")

```julia
sys = init_system(db)
```

```terminal
System:
        Elements: 2
        Phases: 7
```

Point calculation of an equilibrium state:

```julia
X = [0.3, 0.7]
T = 700 # K
equilib(sys, X, T)
```

```terminal
Gamma: 0.02890 mol
        sublattice 1
                Cu: 1.0000
                Zn: 0.0000
        sublattice 2
                Cu: 0.7363
                Zn: 0.2637
        sublattice 3
                Cu: 0.0007
                Zn: 0.9993
Eps: 0.24864 mol
        sublattice 1
                Cu: 0.2269
                Zn: 0.7731
```

## TODO

### Upcoming
- Built-in calculation of equilibrium solidification (temperature dependency of phase constitution)

### Near-future
- Tests
- Documentation

### Faraway-future
- Phase diagram calculation