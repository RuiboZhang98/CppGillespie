# Sequential mutation of cancer cells

Cancer is a genetic disease that results from accumulation of driver mutations which confer a selective growth advantage to tumor cells. To facilitate mathematical quantification of the carcinogenic process, stochastic models can be used to model the accumulation of driver mutations, in particular population sizes and arrival time distributions for premalignant and malignant subpopulations.

# Multi-type branching process

We consider $q$ types of cancer cells. The population vector at given time $t$ is denoted by 

$$\mathbf{N}(t) = (N_0(t), N_1(t),\cdots,N_{q-1}(t)).$$

The type $i$ individuals have birth rate $r_i \geq 0$, death rate $d_i \geq 0$, and $u_{ij}, j\neq i$ which is the transition rate to type $j$. We define the net growth rate as $\lambda_i = r_i - d_i$.

On the level of each individual, we can tell that a type $i$ individual has an exponentially distributed lifespan with parameter $a_i := r_i + d_i + \sum_{j \neq i}u_{ij}$. When a type-$i$ dies, it splits into several new individuals, not necessarily of its own type. We refer to the distribution of newly born individuals as the *offspring distribution* 
$\mathbf{\xi_i} = ( \xi_{i1},\xi_{i2},\cdots,\xi_{iq} )$. 

The probability mass function of $\xi_i$ is

$$
\begin{aligned}
    \mathbb{P}(\mathbf{\xi_i} = \mathbf{0}) &= \frac{d_i}{a_i}, \\
    \mathbb{P}(\mathbf{\xi_i} =  \mathbf{e_i} + \mathbf{e_j}) &= \frac{u_{ij}}{a_i} \quad (j \neq i), \\
    \mathbb{P}(\mathbf{\xi_i} = 2\mathbf{e_i}) &= \frac{b_i}{a_i}.
\end{aligned}
$$

In the above expression, we have set the mutation to be *mutation at division*. This is different from *migration* or *pure mutation* where 

$$\mathbb{P}(\mathbf{\xi_i} = \mathbf{e_j}) = \frac{u_{ij}}{a_i} \quad (j \neq i).$$

# Gillepie Algorithm for Monte Carlo Simulation of multi-type branching processes

Roughly speaking, the algorithm for simulating a single realization of the stochastic process reads:

0. Specify the intial population vector *v0 <- N0*, and the intial time *t <- 0*
1. Multiple each entry of the *v0* by the life span parameter vector $a_i$. Let *a* be the vector of $a_i$.
   *weight <- v0 elementsie product a* 
2. Sample an exponential distribution *Exp* with paramter *weight*. The time for next population change event is
   *t <- t + Exp*
3. Sample an uniform distribution *U*. Use a for loop to figure out the event type *event*. Then update the population vector
   *v0 <- v0 + event*
4. Go back to 1 unless *t > tmax*. 
