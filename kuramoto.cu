#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

#define h 0.01

#define lambda 0.5

#define T 100

#define N 100

__global__ void kuramoto(int *A, double *theta, double *w, double *k, double *prevk, int iter, double adder) {
	int id = threadIdx.x;
	
	// first set it equal to the initial omega values
	k[id] = w[id];
	
	for(int i = 0; i < N; i++) {
		k[id] += lambda * A[i * N + id] * sin(theta[iter * N + i] - (theta[iter * N + id] + (adder * prevk[id])));
	}
}

__global__ void update_theta(double *theta, double *k1, double *k2, double *k3, double *k4, int iter) {
	int id = threadIdx.x;
	
	theta[(iter + 1) * N + id] = theta[iter * N + id] + (h * (k1[id] + 2 * k2[id] + 2 * k3[id] + k4[id]) / 6);
	
	// make sure theta stays between 0 and 2pi
	while(theta[(iter + 1) * N + id] > M_PI) {
		theta[(iter + 1) * N + id] -= 2 * M_PI;
	}
	
	while(theta[(iter + 1) * N + id] < -M_PI) {
		theta[(iter + 1) * N + id] += 2 * M_PI;
	}
}

int main() {

	std::ifstream adj_file("./adjacency.txt");
	std::ifstream phase_file("./initial_phase.txt");
	std::ifstream omega_file("./omega.txt");
	
	const int timesteps = T / h;
	
	double theta[timesteps][N];
	
	double omega[N];
	
	double k0[N];
	
	double *dtheta;
	double *domega;
	double *dk0, *dk1, *dk2, *dk3, *dk4;
	double r[timesteps];
	
	double c = 0, s = 0;

	int A[N][N];
	
	int* dA;
	
	
	// read the contents of adjacency matrix file to A[N][N]
	if(adj_file.is_open()) {
		for(int i = 0; i < N; i++) {
			std::string line;
			std::getline(adj_file, line);
			std::istringstream iss(line);
			for(int j = 0; j < N; j++) {
				iss >> A[i][j];
			}
		}
		adj_file.close();
	}
	
	// read the contents of phase values into theta[N][0]
	if(phase_file.is_open()) {
		int i = 0;
		std::string line;
		while(std::getline(phase_file, line)) {
			std::istringstream iss(line);
			iss >> theta[0][i];
			i++;
		}
		phase_file.close();
	}
	
	// read the contents of omega values into omega[N]
	if(omega_file.is_open()) {
		int i = 0;
		std::string line;
		while(std::getline(omega_file, line)) {
			std::istringstream iss(line);
			iss >> omega[i];
			i++;
		}
		omega_file.close();
	}
	
	// allocate device memory into dA, dtheta, domega
	cudaMalloc(&dA, N * N * sizeof(int));
	cudaMalloc(&dtheta, N * (T / h) * sizeof(double));
	cudaMalloc(&domega, N * sizeof(double));

	
	// allocate device memory for dk0, dk1, dk2, dk3, dk4
	cudaMalloc(&dk0, N * sizeof(double));
	cudaMalloc(&dk1, N * sizeof(double));
	cudaMalloc(&dk2, N * sizeof(double));
	cudaMalloc(&dk3, N * sizeof(double));
	cudaMalloc(&dk4, N * sizeof(double));
	
	// initialize the value of dk0 using k0 from host to device
	for(int i = 0; i < N; i++) {
		k0[i] = 0;
	}
	cudaMemcpy(dk0, k0, N * sizeof(double), cudaMemcpyHostToDevice);
	
	// initialize the values of dtheta using theta from host to device
	cudaMemcpy(dtheta, theta, N * (T / h) * sizeof(double), cudaMemcpyHostToDevice);
	
	// print the values of adjacency matrix as read from the file
	std::cout << "\nAdjacency matrix: \n";
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "\n";
	}
	
	// copy the A values into dA from host to device
	cudaMemcpy(dA, A, N * N * sizeof(int), cudaMemcpyHostToDevice);
	
	// print the values of phase as read from the file
	std::cout << "\nPhase values: ";
	for(int i = 0; i < N; i++) {
		std::cout << theta[0][i] << " ";
	}
	
	// print the values of omega as read from the file
	std::cout << "\nOmega values: ";
	for(int i = 0; i < N; i++) {
		std::cout << omega[i] << " ";
	}
	
	// copy the omega values to domega from host to device
	cudaMemcpy(domega, omega, N * sizeof(double), cudaMemcpyHostToDevice);
	
	std::cout << "\n";
	
	float t = 0;
	int iter = 0;
	
	while(t < T) {
		// run the kuramoto kernel for k1
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk1, dk0, iter, 0);
		cudaDeviceSynchronize();
		
		// run the kuramoto kernel for k2
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk2, dk1, iter, h/2);
		cudaDeviceSynchronize();
		
		// run the kuramoto kernel for k3
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk3, dk2, iter, h/2);
		cudaDeviceSynchronize();
		
		// run the kuramoto kernel for k4
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk4, dk3, iter, h);
		cudaDeviceSynchronize();
		
		// update the value of theta with the values of k1, k2, k3, k4 using the Runge-Kutta formula
		update_theta<<<1, N>>>(dtheta, dk1, dk2, dk3, dk4, iter);
		cudaDeviceSynchronize();
		
		t += h;
		
		iter++;
	}
	
	cudaMemcpy(theta, dtheta, N * timesteps * sizeof(double), cudaMemcpyDeviceToHost);
	
	t = 0;
	iter = 0;
	while(t < T) {
		c = 0;
		s = 0;
		for(int i = 0; i < N; i++) {
			c += cos(theta[iter][i]);
			s += sin(theta[iter][i]);
		}
		r[iter] = sqrt(c*c + s*s) / N;
		t += h;
		iter += 1;
	}
	
	
	// printing the theta values at a certain timestep
	std::cout << "\nTheta values at the 1000th place: \n";
	for(int i = 0; i < N; i++) {
		std::cout << theta[1000][i] << " ";
	}
	
	std::ofstream order_par_file("./order_par.txt");
	
	if(order_par_file.is_open()) {
		for(int i = 0; i < 2000; i += 1) {
        		order_par_file << r[i] << std::endl;
        	}
        }
	
	std::cout << "\n";
	
	cudaFree(dA);
	cudaFree(dtheta);
	cudaFree(domega);
	cudaFree(dk0);
	cudaFree(dk1);
	cudaFree(dk2);
	cudaFree(dk3);
	cudaFree(dk4);
	
	return 0;
}


