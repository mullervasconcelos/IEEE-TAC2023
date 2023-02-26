# Matlab Code for "Robust one-shot estimation over shared networks in the presence of denial-of-service attacks"


# Proactive Jammer

## p2p channel

* **OptimalJammingProbability.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer. 
* **OptimalJammingProbability_VS_d.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer as a function of d. 
* **OptimalJammingProbability_VS_var.m**: Optimal jamming probabilities $\varphi^\star$ as a function of $\sigma^2$. 
* **OptimalJammingProbability_VS_cd.m**: Optimal jamming probability $\varphi^\star$  for the proactive jammer as a function of  $c$ and $d$. 

## Large-scale network
* **OptimalTransmission_Jamming_Policy.m**: Optimal Optimal transmission policy and jamming probability $\varphi^\star$ for the proactive jammer over large-scale network. 
* **OptimalJammingProbability_VS_d.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer over large-scale network as a function of d. 
* **OptimalJammingProbability_VS_kappa.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer over large-scale network as a function of $\bar{\kappa}$. 
* **OptimalJammingProbability_VS_cd.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer over large-scale network as a function of  c and d. Here, X ~ N(0,1). 
* **OptimalJammingProbability_VS_ckappa.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer over large-scale network as a function of  c and $\bar{\kappa}$. Here, X ~ N(0,1). 
* **OptimalJammingProbability_VS_dkappa.m**: Optimal jamming probability $\varphi^\star$ for the proactive jammer over large-scale network as a function of d and $\bar{\kappa}$. Here, X ~ N(0,1). 



# Reactive Jammer

## 1 Dimension Gaussian variable

### functions
* **MonteCarlo_AlgorithmComparison.m**: Monte Carlo simulations for PGA-CCP and GDA in 1 dimension
* **Optimalalphabeta_VS_var**: Optimal jamming probabilities $\alpha^\star$ and $\beta^\star$ as a function of $\sigma^2$. 
* **PGA_CCP_Algorithm**: Convergence of PGA-CCP 
* **GDA_Algorithm.m**: Convergence of GDA 

### auxiliary functions

* **FirstNashEqulibiumChecker.m**: Check whether approximate First Nash Equilibrium is satisfied
* **grad_PGA.m** Gradients for PGA
* **grad_CCP.m**: Gradients for CCP
* **grad_GD.m**: Gradients for GD


## Multidimensional Guassian vector
**Description:** the algorithm is accelerated by replacing iterations with matrix computation.
### functions
* **MonteCarlo_MultiDim_AlgorithmComparison.m**: Monte Carlo simulations for PGA-CCP and GDA for multidimensional state, the expectation is estimated by using the average of random samples
* **MultiDim_PGA_CCP_Algorithm** Convergence of PGA-CCP for multidimensional state
* **MultiDim_GDA_Algorithm.m**: Convergence of GDA for multidimensional state

### auxiliary functions
* **batch_grad_PGA.m**  Estimation of Gradients for PGA
* **batch_grad_CCP.m**: Estimation of Gradients for CCP
* **batch_grad_GD.m**:  Estimation of Gradients for GD
* **batch_FirstNashEqulibiumChecker.m**: Check whether the estimated approximate First Nash Equilibrium is satisfied