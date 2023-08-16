# **Simulating the 3D Ising Model Using Markov-Chain Monte Carlo**

# Introduction

**Goal:**\
To compute the transition/critical/Curie temperature for the model, which is the temperature at which the system exhibits a (second order) phase transition from a ferromagnetic phase to a paramagnetic phase.

**Procedure:**
1. simulate the 3-dimensional Ising Model using the Metropolis-Hasting's algorithm at various physical temperatures.
2. compute physical observables for each temperature.
3. observe the behaviour as a function of temperature. 

This yields the value of the critical temperature, which can be seen as the temperature corresponding to a rapid increase/decrease/discontinuity in the plots of the physical observables.


# Theoretical Background

## A Brief Review of the Physics of Phase Transitions

**Phase:**
- a region of parameters space where the properties of the system are analytic, Taylor sums converge, perturbation theory works well, the partition function is well defined etc.

**Phase boundaries:**
- parameter values where the aforementioned properties are not smooth and adopt discontinuities.
- the large particle number limit of thermodynamics breaks down at phase boundaries and fluctuations build up.

**Phase transition:**
- an abrupt, discontinuous change in the properties of a system. 
also associated with a symmetry breaking.
- loosely speaking, a phase transition is a transition between ordered and disordered regimes. For example ice gains rotational symmetry on transitioning to water

According to the conventional **classification of phase transitions**, a transition is
- **first-order** if first order quantities (for example energy and magnetization) w.r.t the partition function are discontinuous with respect to the system's order parameter (in this case, the temperature). 
- **second-order** if first order quantities are continuous, but second order quantities (for example heat capacity and susceptibility) w.r.t the order parameter are discontinuous.

**Critical points:**
- points where phase boundaries vanish. These points are characterized by attributes like universality and scale invariance.


## Ising Model

1. A mathematical model of binary states on a lattice, that does not correspond to any one particular  physical system.
    - for $N$ points per dimension and $D$ dimensional model, the # of possible configurations is $2^{N^D}$ = $10^{N^D \log2 }$ 
    - can represent a magnetic system with electron spins, beta-brass, a lattice gas etc.

2. The Hamiltonian (energy) of the system is given by the following expresion,
$$\mathcal{H} = -J\sum_{\text{lattice }(i)}\sum_{\text{NN} (j)} s_i s_j - B \sum_{\text{lattice }(i)}s_i \tag{1}$$

where $B \sim$ magnetic field and $J \sim$ spin interaction term. The index $i$ runs over all lattice points and $j$ runs over the nearest neighbors of a lattice point $i$ (every point in a 3D lattice has $6$ nearest neighbors).

Ferromagnetism arises when a collection of atomic spins align such that their associated magnetic moments all point in the same direction, yielding a net magnetic moment which is macroscopic in size. The simplest theoretical description of ferromagnetism is given by Ising model.

4. The model was given by Ernst Ising in 1925, who also completely solved the 1D version of the model. The 2D model was completely solved by Onsager (for $B=0$), uisng transfer matrices, and is well known for being extremely complicated. There are no analytical solutions for higher dimensions, which is what motivates the need for a computational method.


## Metropolis Hasting's Algorithm

This algorithm generates a markov chain, where each element in the chain is a state of the lattice. The following are the necessary properties for this algorithm to work,
1. *Time homogeniety:* the transition probability of state $i$ to $i+1$ is independent of state of the markov chain.
2. *Stationarity of final distribution:* application of the transition probability matrix to the stationary distribution vector returns the same vector, i.e its an eigenvector with eigenvalue $1$.
3. *Irreducibility of transition probabilities:* There exists a time after which any possible state in the state space is arrived at from any other state in the space with a non-zero probability.
4. *Aperiodicity of the chain:* the greatest common divisor of all times after which the chain returns a state to itself is $1$.

If these properties hold, then monte carlo sums converges to the correct expectation values.

**Physical Setup:**
- The lattice is placed in a (much larger) reservoir at temperature  T . This ensures that the probability distribution followed by the system is the Boltzmannian distribution (  x∼exp(−x) )
- Periodic boundary conditions are imposed, which means that the lattice is a 3-torus

**Implementation of the MH Algo:** To build a chain of lattice configurations $\{S_i\}$.
From lattice state $S_i$, $S_{i+1}$ is arrived at using the following algorithm,
1. Pick one lattice point at random.
2. Calculate the change in energy $\Delta E$ for swapping the spin at that point.
3. If $\Delta E < 0$, then flip spin and the hence formed lattice is the state $S_{i+1}$.
4.  If $\Delta E > 0$, then flip spin with probability $exp(-\beta \Delta E)$, and the hence formed lattice is the state $S_{i+1}$, where $1/\beta = kB \cdot T$.


## Observables Extracted from the Simulation

The following will be evaluated using our simulations ($N_{MC}$ = total Monte Carlo time steps, $N_{lat}$ = lattice dimensions):

1. **Average Energy per lattice point**: $$\frac{<\mathcal{H}>}{N_{lat}^3} = \frac{1}{N_{lat}^3 \times N_{MC}} \sum_{\text{MC}(t)} \mathcal{H}_t \tag{2}$$

2. **Average Magnetization per lattice point**:

$$\frac{<\mathcal{M}>}{N_{lat}^3} = \frac{1}{N_{lat}^3 \times N_{MC}} \sum_{\text{MC}(t)} \mathcal{M}_t \tag{3}$$

$$\mathcal{M} = \sum_{\text{lattice }(i)}s_i \tag{4}$$

3. **Heat Capacity per lattice point**: 
$$\frac{<\mathcal{C}>}{N_{lat}^3} = \frac{1}{N_{lat}^3 \times N_{MC}} \sum_{\text{MC}(t)} \mathcal{C}_t \tag{5}$$ where $$\mathcal{C} = k_B\beta^2 \left(<\mathcal{H}^2> - <\mathcal{H}>^2\right) \tag{6}$$

4. **Succeptibility per lattice point**: 
$$\frac{<\mathcal{\chi}>}{N_{lat}^3} = \frac{1}{N_{lat}^3 \times N_{MC}} \sum_{\text{MC}(t)} \mathcal{\chi}_t \tag{7}$$ where $$\mathcal{\chi} = \beta\left(<\mathcal{M}^2> - <\mathcal{M}>^2 \right) \tag{8}$$

**Note:** These quantities are actually evaluated over $N_{MC} - N_{therm}$ times at a gap of $N_{corr}$ states from the chain ($N_{corr}$ = correlation time steps, $N_{therm}$ = thermalizaion time steps). This is done to ensure that only sum those states in the chain after the chain has thermalized, and only over uncorrelated states in the chain.


# Implementation

## Simulation Parameters

1. $B = 0$ $\equiv$ magnetic field * |magnetic moment|

    The magnetic field is set to 0. This is done to observe for spontaneous magnetization even in the absence of magnetic field.

2. $J = 1 \equiv$ interaction energy parameter * |magnetic moment|**2 and $kB = 1$.

    This implies that temperature and energy both are dimensionless, and hence the value of the transition temperature evaluated is units invariant.

3. $N_{lat}$ = {15,20}

    The simulation was run for latttice sizes 15 and 20

4. $N_{chain}$ = 1e6

    Number of MC steps in each chain.

5. $N_{proc}$ = 10

    Number of chains for each temperature.
    
    
## Flow of Information in Scripts

1. A set of config files are passed into the master script **run_ising_simulation.py** as input. These contain information about the simulation parameters.

2. In the master script, $N_{m-procs}$ config files lead to a launch of a master process in parallel. 

3. Each master process then generates $N_{procs}$ children chains parallely as children processes for temperatures in np.linspace(T_start, T_end, N_T) serially.

3. The master scripts uses functions in ***comp_quan.py***, ***setup_lat.py*** and ***sim_ising.py***, which were all cythonized to acheive a speedup.

4. The output of this is a set of hdf files for each chain for each temperature for each config file in a correspondingly located folder in the subdirectories of folder containing the master script. This output file contains the value of the energy (hamiltonian) and magnetization of the lattice at each MC step.

5. The output files for chains corresponding to one temperature are averaged, and the output of this is passed into *comp_ising_quant.py*. This computes the observables listed in the previous section. These observables are saved in a json file.

6. The *plot_ising_quant.py* script reads data from this json file and plots the data. These plots are present in the following section.


The following are flow-charts that summarize the flow of information in the simulation and computation scripts.
![image](https://github.com/a-ravichandran/3d_ising_model/assets/90862509/f5bd394f-a5e9-47f5-86cc-cfa3cd7905fa)
![image](https://github.com/a-ravichandran/3d_ising_model/assets/90862509/0bc6a1dd-338c-4f2b-bf5e-8dbb24762ffd)


# Observations, Results and Conclusions

## MC History
The two plots below indicate the folowing:

1. All chains thermalize eventually, but the chains corresponding to low temperatures take longer to thermalize.

2. Thermalization was also observed to be faster in the smaller lattice. This is because the state space of a larger lattice is larger.

3. Hence the N_therm, and hence the number of initial MC steps discarded and not included in the computation of observables was separately set for different temperatures and the two lattices.

![image](https://github.com/a-ravichandran/3d_ising_model/assets/90862509/8618a88b-9205-486b-90f3-38c139ced35f)
![image](https://github.com/a-ravichandran/3d_ising_model/assets/90862509/62a20a39-eacd-4011-b5c6-9d28ee7d3802)



## Observables as a function of temperature

1. The energy and magnetization show a continuous behavior in the entire temperature range, with a continuous jump observable in around $T \sim 4.5$. 

2. At low temperatures, the average magnetization per cell is close to $1$, which indicates that most spins on the lattice are aligned to each other. Which implies that the total magnetization of the system is large i.e. this phase is ferromagnetic. This is also characterized by more negative energy since the aligned spins reduce the hamiltonian value as per equation (1).

3. On moving to higher temperatures, the average magnetization drops to $0$ suddenly, which indicates that the spins on the lattice randomly cancel each other out, which implies that the total magnetization is very low. This is the paramagnetic phase. This is also characterized by a larger hamiltonian value from equation (1). 

4. By studying the heat capacity and succeptibility plots from higher to lower temperature, show a jump around $T \sim 4.5$, which is the critical temperature for the 3D Ising model for $J=1$ and $kB=1$.

5. The aforementioned jump is higher for the larger lattice, which is the expected behavior, since the extent of the jump scales with the size of the lattice.

6. Since it is the second order quantities that are discontinuous, the loss of spontaneous magnetization in a ferromagnetic material as the temperature exceeds the critical temperature is a second-order phase transition.

![image](https://github.com/a-ravichandran/3d_ising_model/assets/90862509/42e470e7-4eef-4a45-b8a8-443f84bde496)


# Limitations and Possible Future Work
**Low acceptance probability:**
- The markov-chains at low temperatures takes a larger number of steps to thermalize compared to the ones at higher temperatures.

**Phase transitions only occur in infinite systems:**
- Because delta functions usually scale with the spacial volume of a system. But indications of their presence can be seen in finite systems, which we did see in the plots.

**Algorithm is slow:**
- There are a few faster alternatives for example the **heat bath algorithm**, the **Wolff cluster** algorithm.

**Used a coarse temperature mesh:**
- Since one knows now that it is  ∼4.5 , one could have a finer grid of temperature around that value to better estimate it. 
- One can also study the critical behaviour around that temperature, to extract the critical exponents like in the Curie-Weiss Law for magnetic susceptibility.

**B=0:**
- I set  B=0  to study spontaneous magnetization, but the same set of scripts can be used to simulate the model at  B≠0  to study things like Hysterisis, where the behavior of the system depends on the history of the magnetic field.
