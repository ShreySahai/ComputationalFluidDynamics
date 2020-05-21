# Need for Upwind Scheme: A visualisation

Need for upwind scheme in numerical solutions for advection-diffusion type problems is illustrated with the solution of 1-D steady state axial convection-diffusion slug flow problem.The below library shows the exact solution and the Finite Difference Scheme solution for different ranges of peclet number


## Getting started 

The repository can be cloned and the library `Visualise.py` can be simply executed for the visualisation

### Prerequisites 

The library uses the following two python libraries

* numpy
* math
* matplotlib

### Problem statement

One Dimensional steady state diffusion-convection flow with a temperature different along some length can be described by the energy equation as 

![Governing equation and supporting diagram](/NeedForUpwindScheme_Visualisation//images/1.png)

Non dimensionalising and Central difference scheme leads to Tridiagonal formulation

(0.5+0.25Pec)Ti-1 -Ti +(0.5-0.25Pec)Ti+1=0 for i=1 to N-1

Ti=0=T0 and Ti=N=T1

where as the exact solution is 

T=(exp(Pex)-1)/(exp(p)-1)

where Pec =Pe△x 

The interval x=(0,10) is broken with 100 nodes. The resulting temperature profile is normalised.


## Usage

To be updated soon

![The solution for different Peclet number](/NeedForUpwindScheme_Visualisation//images/Example.png)


## Observations and Discussion

* The distribution is linear at Pec<<1. Signifying insignificance of advection term
* The temperature difference concentrates to a thin region as Pe increases signifying increasing effect on advection. The advection term tends to reduce the temperature variance where as the diffusion here is acting as a source term and acts in to increase the temperature difference across the length
* Oscillations  oscillations in the numerical solution starts appearing at high Peclet number (Pec>=2) signifying the error in formulations and suggesting upwind scheme as a solution


