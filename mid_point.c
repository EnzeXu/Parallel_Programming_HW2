#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.14159265358979323846 // accurate pi
#define N_MAX 1100000 // maximum of N
#define X_MIN 0.0
#define X_MAX 1.0

double func(double x) {
    return 4.0 / (1.0 + x * x);
}

// double originFunc(double x) {
//     return 4.0 * atan(x);
// }

// double truth(void) {
//     return originFunc(X_MAX) - originFunc(X_MIN);
// }

double sequentialTest(int n) {
    int i;
	double sum = 0;
	for (i = 0; i < n; i++) {
        double tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
        sum += tmp;
	}
	return sum;
}

double parallelTestCritical(int n) {
	int i;
    double sum = 0;
	// omp_set_num_threads(n_thread);
	#pragma omp parallel for schedule(static)
	for (i = 0; i < n; i++) {
        double tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
        #pragma omp critical
        sum += tmp;
	}
	return sum;
}

double parallelTestAtomic(int n) {
	int i;
    double sum = 0;
	// omp_set_num_threads(n_thread);
	#pragma omp parallel for schedule(static)
	for (i = 0; i < n; i++) {
        double tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
        #pragma omp atomic
        sum += tmp;
	}
	return sum;
}

double parallelTestReduce(int n) {
	int i;
    double sum = 0;
	// omp_set_num_threads(n_thread);
	#pragma omp parallel for schedule(static) reduction(+:sum)
	for (i = 0; i < n; i++) {
        double tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
        sum += tmp;
	}
	return sum;
}


void driver(void) {
	// plot 1.1 - Power
	int nListPow[10] = {2, 4, 8, 16, 32, 64, 128, 256, 51200, 1024};
	printf("Plot 1.1 - Power\n");
	printf("      n    seq error      par error\n");
	for (int i = 0; i < 10; i++){
		double seqError = fabs(sequentialTest(nListPow[i]) - PI) / PI;
		double parError = fabs(parallelTestReduce(nListPow[i], 50) - PI) / PI;
		printf("%7d    %.9lf    %.9lf\n", nListPow[i], seqError, parError);
	}
	printf("\n");

	// plot 1.2 - Linear
	int nListLinear[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
	printf("Plot 1.2 - Linear\n");
	printf("      n    seq error      par error\n");
	for (int i = 0; i < 10; i++){
		double seqError = fabs(sequentialTest(nListLinear[i]) - PI) / PI;
		double parError = fabs(parallelTestReduce(nListLinear[i], 50) - PI) / PI;
		printf("%7d    %.9lf    %.9lf\n", nListLinear[i], seqError, parError);
	}
	printf("\n");


    int n = 100000000;
    double tSeqStart = clock();
	sequentialTest(n);
	double tSeq = clock() - tSeqStart;
    printf("seq time baseline = %lf s\n", tSeq / CLOCKS_PER_SEC);

	// plot 2.1
	printf("Plot 2.1 - Critical\n");
	printf("n_thread  efficiency\n");
	for (int i = 5; i <= 50; i+=5){
        omp_set_num_threads(i);
		double t0 = clock();
		parallelTestCritical(n);
		double t1 = clock();
		double efficiency = tSeq / (t1 - t0);
		printf("%6d    %.15lf\n", i, efficiency);
	}

    // plot 2.2
	printf("Plot 2.2 - Atomic\n");
	printf("n_thread  efficiency\n");
	for (int i = 5; i <= 50; i+=5){
        omp_set_num_threads(i);
		double t0 = clock();
		parallelTestAtomic(n);
		double t1 = clock();
		double efficiency = tSeq / (t1 - t0);
		printf("%6d    %.15lf\n", i, efficiency);
	}

    // plot 2.3
	printf("Plot 2.3 - Reduce\n");
	printf("n_thread  efficiency\n");
	for (int i = 5; i <= 50; i+=5){
        omp_set_num_threads(i);
		double t0 = clock();
		parallelTestReduce(n);
		double t1 = clock();
		double efficiency = tSeq / (t1 - t0);
		printf("%6d    %.15lf\n", i, efficiency);
	}
	return;
}

int main(int argc, char **argv) {
    // printf("critical: %.15lf\n", parallelTestCritical(1000000, 10));
    // printf("atomic: %.15lf\n", parallelTestAtomic(1000000, 10));
    // printf("reduce: %.15lf\n", parallelTestReduce(1000000, 10));
    driver();
    return 0;
}