/*
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
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>

#define h 0.01

#define lambda 0.1

__global__ void kuramoto(double *A, double *theta, double *w, 
				double *k, double *prevk, int N, int iter, double adder) {
	int id = threadIdx.x;
	
	theta[id] = w[id];
	
	for(int i = 0; i < N; i++) {
		k[id] += lambda * A[i*N + id] * sin((theta[i] + adder * prevk[i]) - theta[id + iter*N]);
	}
}

__global__ void update_theta(double *theta, double *k1, double *k2, double *k3, double *k4, int iter, int N) {
	theta[(iter+1)*N + threadIdx.x] = theta[iter*N + threadIdx.x] + (double)(h * (k1 + 2*k2 + 2*k3 + k4)) / 6;
}

int main() {
	FILE *fptr_1, *fptr_2, *fptr_3;

	int N = 0;
	
	fptr_1 = fopen("./adjacency.txt", "r");
	fptr_2 = fopen("./initial_phase.txt", "r");
	fptr_3 = fopen("./omega.txt", "r");
	
	double ph, ff;
	
	unsigned float t = 0;
	
	do {
		ff = fscanf(fptr_2, "%lf", &ph);
		N++;
	} while(ff != EOF)
	
	double *theta, *dtheta;
	theta = (double *)malloc(N * 10000 * sizeof(double));
	cudaMalloc(&dtheta, N * 10000 * sizeof(double));
	
	for(int i = 0; i < N; i++) {
		fscanf(fptr_2, "%lf", theta[i]);
	}
	
	cudaMemcpy(dtheta, theta, N * 10000 * sizeof(double), cudaMemcpyHostToDevice);
	
	int *A, *dA;
	A = (int *)malloc(N * N * sizeof(int));
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			fscanf(fptr_1, "%d", A[i*N + j]);
		}
	}
	
	cudaMemcpy(dA, A, N * N * sizeof(int), cudaMemcpyHostToDevice);
	
	double *omega, *domega;
	
	for(int i = 0; i < N; i++) {
		fscanf(fptr_3, "%lf", omega[i]);
	}
	
	cudaMemcpy(domega, omega, N * sizeof(double));
	
	int iter = 0;
	
	double *dk1, *dk2, *dk3, *dk4;
	cudaMalloc(&dk1, N*sizeof(double));
	cudaMalloc(&dk2, N*sizeof(double));
	cudaMalloc(&dk3, N*sizeof(double));
	cudaMalloc(&dk4, N*sizeof(double));
	
	while(t < 100) {
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk1, NULL, N, iter, 0);
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk2, dk1, N, iter, h/2);
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk3, dk2, N, iter, h/2);
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk4, dk3, N, iter, h);
		
		update_theta<<<1, N>>>(dtheta, dk1, dk2, dk3, dk4, iter, N);
		
		t += h;
		
		iter++;
	}
	
	fclose(fptr_1);
	fclose(fptr_2);
	fclose(fptr_3);
	
	return 0;

}


