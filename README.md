kuramoto oscillators

equation:
	
	d(theta[i])/dt = w[i] + lambda * sum( A[i][j] * sin(theta[j] - theta[i])
	
	solving the above using rk4 method:
	
	let df[i] = d(theta[i])/dt = f(theta[i]);
	
	then k1 = f(theta[i]);
	
	and k2 = f(theta[i] + h*k1/2)
	
	k3 = f(theta[i] + h*k2/2)
	
	k4 = f(theta[i] + h*k3)
	
	then we can find theta[i + 1] using the following formula:
	
	theta[i+1] = theta[i] + h*(k1 + 2*k2 + 2*k3 + k4) / 6

