# Multi-type branching process

We consider $q$ types of individuals. The population of individuals is denoted by 

$$\mathbf{N}(t) = (N_0(t), N_1(t),\cdots,N_{q-1}(t)).$$

The type-$i$ individuals have birth rate $r_i \geq 0$, death rate $d_i \geq 0$, and $u_{ij}, j\neq i$ which is the transition rate to type-$j$. We define the net growth rate as $\lambda_i = r_i - d_i$.

On the level of each individual, we can tell that a type-$i$ individual has an exponentially distributed lifespan with parameter $a_i := r_i + d_i + \sum_{j \neq i}u_{ij}$. When a type-$i$ dies, it splits into several new individuals, not necessarily of its own type. We refer to the distribution of newly born individuals as the *offspring distribution* 
$\mathbf{\xi}_i= (\xi_{i1},\xi_{i2},\cdots,\xi_{iq})$. The probability mass function of $\xi_i$ is

$$
\begin{aligned}
    \mathbb{P}(\mathbf{\xi}_{i} = \mathbf{0}) &= \frac{d_i}{a_i}, \\
    \mathbb{P}(\mathbf{\xi}_{i} =  \mathbf{e}_i + \mathbf{e}_j) &= \frac{u_{ij}}{a_i} \quad (j \neq i), \\
    \mathbb{P}(\mathbf{\xi}_{i} = 2\mathbf{e}_i) &= \frac{b_i}{a_i}.
\end{aligned}
$$

In the above expression, we have set the mutation to be *mutation at division*. This is different from *migration* or *pure mutation* where 

$$\mathbb{P}(\mathbf{\xi}_{i} = \mathbf{e}_j) = \frac{u_{ij}}{a_i} \quad (j \neq i).$$
