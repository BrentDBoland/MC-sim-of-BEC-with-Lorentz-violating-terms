# MC-sim-of-BEC-with-Lorentz-violating-terms
A java-based Monte Carlo algorithm to solve for the ground state of a Bose-Einstein condensate of rubidium atoms trapped in a pseudo 1D potential modeled by the Gross-Pitaevskii equation with Lorentz Violating terms from the Standard Model Extension added to the Hamiltonian. Built using opensourcephysics package. The program uses random walkers to estimate the energy expectation value, <E> of a wave function constructed from a linear combination of eigen functions of the harmonic oscillator. After determining <E>, the program iteratively perturbs the basis vector (coefficients of the eigen functions) of the wave function and accepts changes that reduce <E>. GPEwithLVtermBECgroundstateMonteCarloApp.jar repeats this process for an iteratively increasing strength in a chosen Lorentz violating term.
# Executables
Two jar files included. HarmOscBECgroundstateMonteCarloApp.jar solves for the wave function of with minimum energy state for the simple harmonic oscillator. GPEwithLVtermBECgroundstateMonteCarloApp.jar solves for the wave function of with minimum energy state for a simple harmonic oscillator with a Gross-Pitaevskii interaction term and two Lorentz violating terms.
# GUI
<img width="353" height="557" alt="image" src="https://github.com/user-attachments/assets/a0eb3835-a606-495d-b8ca-763e70d88dcb" />

The GUI allows the user to set 
1) The mass of the atoms.
2) The strength of the trap omega.
3) The boson–boson s-wave scattering length.
4) The number of harmonic oscillator eigen states to use (limits the dimensions of the basis vector dimensions to the first n lowest energy states).
5) The number of samples to use for the <E> computation.
6) The width of domain in space to consider (the farther from x=0 the less the contribution to <E> so limiting the domain can reduce computation time).
7) The norm of the wavefunction (in reality this should always be 1, but scaling it can help avoid rounding errors).
8) The maximum acceptance rate to terminate the search.
9) The norm of the step size for perturbations in the basis vector.
10) The step size for the random walkers used in calculating <E>.
11) the Lorentz violating term to sweep through (0 for A and 1 for F).
12) Set the value of the other Lorentz violating term (Ex: if the above number is 0 and this number is 0 then F=0).
13) Allows the user to choose an initial test function that is one of the n basis vectors.
<img width="805" height="277" alt="image" src="https://github.com/user-attachments/assets/5b89fb0f-fda5-42b1-ae93-4a0d2e1613bc" />

<img width="806" height="277" alt="image" src="https://github.com/user-attachments/assets/c45a657c-4d83-41dd-b805-669639b20eaf" />

Plots of both the wave function and the Hamiltonian of the wave function are displayed. For an exact ground state the Hamiltonian should be the wave function scaled by E_0.


