#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cuda.h>

#define h 0.01

#define lambda 0.8

__global__ void kuramoto(int *A, double *theta, double *w, double *k, double *prevk, int N, int iter, double adder) {
	int id = threadIdx.x;
	
	k[id] = w[id];
	
	for(int i = 0; i < N; i++) {
		k[id] += lambda * A[i*N + id] * sin(theta[iter*N + i] - (theta[iter*N + id] + (adder * prevk[i])));
	}
}

__global__ void update_theta(double *theta, double *k1, double *k2, double *k3, double *k4, int iter, int N) {
	int id = threadIdx.x;
	theta[(iter+1)*N + id] = theta[iter*N + id] + (double)(h * (k1[id] + 2*k2[id] + 2*k3[id] + k4[id])) / 6;
				
	while(theta[(iter+1)*N + id] > M_PI) {
		theta[(iter+1)*N + id] = 2 * M_PI - theta[(iter+1)*N + id];
	}
	
	while(theta[(iter+1)*N + id] < -M_PI) {
		theta[(iter+1)*N + id] += 2*M_PI;
	}
}

int main() {
	FILE *fptr_1, *fptr_2, *fptr_3;
	
	fptr_1 = fopen("./adjacency.txt", "r");
	fptr_2 = fopen("./initial_phase.txt", "r");
	fptr_3 = fopen("./omega.txt", "r");
	
	double *theta, *dtheta;
	double *omega, *domega;
	double *k0, *dk0, *dk1, *dk2, *dk3, *dk4;
	
	int N = 10;

	int *A, *dA;
	
	A = (int *)malloc(N * N * sizeof(int));
	cudaMalloc(&dA, N * N * sizeof(int));
	
	theta = (double *)malloc(N * 10000 * sizeof(double));
	cudaMalloc(&dtheta, N * 10000 * sizeof(double));
	
	omega = (double*)malloc(N * sizeof(double));
	cudaMalloc(&domega, N * sizeof(double));
	
	cudaMalloc(&dk0, N*sizeof(double));
	cudaMalloc(&dk1, N*sizeof(double));
	cudaMalloc(&dk2, N*sizeof(double));
	cudaMalloc(&dk3, N*sizeof(double));
	cudaMalloc(&dk4, N*sizeof(double));
	
	k0 = (double *)malloc(N * sizeof(double));
	for(int i = 0; i < N; i++) k0[i] = 0;
	cudaMemcpy(dk0, k0, N * sizeof(double), cudaMemcpyHostToDevice);
	
	for(int i = 0; i < N; i++) {
		fscanf(fptr_2, "%lf", &theta[i]);
		printf("%lf ", theta[i]);
	}
	cudaMemcpy(dtheta, theta, N * 10000 * sizeof(double), cudaMemcpyHostToDevice);
	printf("\n");
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			fscanf(fptr_1, "%d", &A[i*N + j]);
			printf("%d ", A[i*N + j]);
		}
		printf("\n");
	}
	cudaMemcpy(dA, A, N * N * sizeof(int), cudaMemcpyHostToDevice);
	
	for(int i = 0; i < N; i++) {
		fscanf(fptr_3, "%lf", &omega[i]);
		printf("%lf ", omega[i]);
	}
	cudaMemcpy(domega, omega, N * sizeof(double), cudaMemcpyHostToDevice);
	
	printf("\n");
	float t = 0;
	int iter = 0;
	
	while(t < 100) {
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk1, dk0, N, iter, 0);
		cudaDeviceSynchronize();
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk2, dk1, N, iter, h/2);
		cudaDeviceSynchronize();
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk3, dk2, N, iter, h/2);
		cudaDeviceSynchronize();
		
		kuramoto<<<1, N>>>(dA, dtheta, domega, dk4, dk3, N, iter, h);
		cudaDeviceSynchronize();
		
		update_theta<<<1, N>>>(dtheta, dk1, dk2, dk3, dk4, iter, N);
		cudaDeviceSynchronize();
		
		t += h;
		
		iter++;
	}
	
	cudaMemcpy(theta, dtheta, N * 10000 * sizeof(double), cudaMemcpyDeviceToHost);
	
	for(int i = 0; i < N; i++) {
		printf("%lf ", theta[9900 + i]);
	}
	
	printf("\n");
	
	fclose(fptr_1);
	fclose(fptr_2);
	fclose(fptr_3);
	
	return 0;
}


