# kuramoto oscillators

### Equation:

```math	
	\frac{d(\theta_i)}{dt} = w_i + lambda*\sum_{j = 1}^{N} A_{ij} \sin(\theta_j - \theta_i)
```
	
Solving the above using rk4 method:

Let `df_i = d(\theta_i)/dt = f(\theta_i)`
	
Then,

```math
	k_1 = f(\theta_i)\newline
	
	k_2 = f(\theta_i + h*\frac{k_1}{2})\n
	
	k_3 = f(\theta_i + h*\frac{k_2}{2})\n
	
	k_4 = f(\theta_i + h*k_3)
```
	
Then we can find `\theta_{i + 1}` using the following formula: `\theta_{i+1} = \theta_i + h*(k_1 + 2*k_2 + 2*k_3 + k_4) / 6`

