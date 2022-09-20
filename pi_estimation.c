#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define RANDOM_N 1000000 // used to generate random number in [0, 1]
#define PI 3.14159265358979323846 // accurate pi
#define N_MAX 1100000 // maximum of N


double getRandom(void) {
	return ((rand() % (RANDOM_N + 1)) / (double) RANDOM_N);
}

int insideQuaterCircle(double x, double y) {
	return (int) (x * x + y * y <= 1);
}

double sequentialTest(int n) {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		sum += insideQuaterCircle(getRandom(), getRandom());
	}
	double pi = sum / (double) n * 4.0;
	// printf("[Sequential] The estimated pi = %lf\n", pi);
	return pi;
}

double parallelTest(int n, int n_thread) {
	// double t0 = clock();
	int record[N_MAX];
	int i;
	// double t1 = clock();
	omp_set_num_threads(n_thread);
	#pragma omp parallel for
	for (i = 0; i < n; i++) {
		record[i] = insideQuaterCircle(getRandom(), getRandom());
	}
	// double t2 = clock();
	int sum = 0;
	for (int i = 0; i < n; i++) {
		sum += record[i];
	}
	double pi = sum / (double) n * 4.0;
	// double t3 = clock();
	// printf("t0-1: %lf t1-2: %lf t2-3: %lf\n", (t1 - t0)/CLOCKS_PER_SEC, (t2 - t1)/CLOCKS_PER_SEC, (t3 - t2)/CLOCKS_PER_SEC);
	// printf("[Parallel] The estimated pi = %lf\n", pi);
	return pi;
}

void driver(void) {
	// plot 1.1 - Power
	int nListPow[20] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};
	// double seqErrorList[20] = {};
	// double parErrorList[20] = {};
	printf("Plot 1.1 - Power\n");
	printf("      n    seq error   par error\n");
	for (int i = 0; i < 20; i++){
		double seqError = fabs(sequentialTest(nListPow[i]) - PI) / PI;
		double parError = fabs(parallelTest(nListPow[i], 44) - PI) / PI;
		// seqErrorList[i] += seqError;
		// parErrorList[i] += parError;
		printf("%7d    %.6lf    %.6lf\n", nListPow[i], seqError, parError);
	}
	printf("\n");

	// plot 1.2 - Linear
	int nListLinear[20] = {50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 550000, 600000, 650000, 700000, 750000, 800000, 850000, 900000, 950000, 1000000};
	printf("Plot 1.2 - Linear\n");
	printf("      n    seq error   par error\n");
	for (int i = 0; i < 20; i++){
		double seqError = fabs(sequentialTest(nListLinear[i]) - PI) / PI;
		double parError = fabs(parallelTest(nListLinear[i], 44) - PI) / PI;
		printf("%7d    %.6lf    %.6lf\n", nListLinear[i], seqError, parError);
	}
	printf("\n");

	// // plot 2
	// int n = 1000000;
	// double tSeqStart = clock();
	// sequentialTest(n);
	// double tSeq = clock() - tSeqStart;
	// printf("Plot 2\n");
	// printf("seq time baseline = %lf s\n", tSeq / CLOCKS_PER_SEC);
	// printf("n_thread  efficiency\n");
	// for (int i = 1; i <= 50; i++){
	// 	double t0 = clock();
	// 	parallelTest(n, i);
	// 	double t1 = clock();
	// 	double efficiency = tSeq / (t1 - t0);
	// 	printf("%6d    %.6lf\n", i, efficiency);
	// }
	return;
}



int main() {
	srand(time(0));
	driver();
	// sequentialTest(1000000);
	// parallelTest(100000, 10);

	return 0;
}
