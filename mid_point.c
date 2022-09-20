#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.14159265358979323846 // accurate pi
#define X_MIN 0.0
#define X_MAX 1.0
#define P_DEFAULT 10

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

double parallelTestCritical(int n, int n_thread) {
	int i;
    double sum = 0;
    double tmp;
    #pragma omp parallel num_threads(n_thread) private(tmp)
    {
        int tid = omp_get_thread_num();
        double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            part_sum += tmp;
        }
        #pragma omp critical
        sum += part_sum;
    }
    return sum;
}

double parallelTestAtomic(int n, int n_thread) {
	int i;
    double sum = 0;
    double tmp;
    #pragma omp parallel num_threads(n_thread) private(tmp)
    {
        int tid = omp_get_thread_num();
        double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            part_sum += tmp;
        }
        #pragma omp atomic
        sum += part_sum;
    }
    return sum;
}

double parallelTestReduce(int n, int n_thread) {
	int i;
    double sum = 0;
    double tmp;
    #pragma omp parallel num_threads(n_thread) private(tmp) reduction(+:sum)
    {
        int tid = omp_get_thread_num();
        // double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            // part_sum += tmp;
            sum += tmp;
        }
        // #pragma omp atomic
        // sum += part_sum;
    }
    return sum;
}


void driver(int n, int pmax) {
	// plot 1.1 - Power
	int nListPow[20] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};
	printf("Plot 1.1 - Power\n");
	printf("n       pi_sequential   pi_parallel     error_seq       error_par\n");
	for (int i = 0; i < 20; i++){
        double seq_res = sequentialTest(nListPow[i]);
		double seq_error = fabs(seq_res - PI) / PI;
        double par_res = parallelTestReduce(nListPow[i], P_DEFAULT);
		double par_error = fabs(par_res - PI) / PI;
		printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", nListPow[i], seq_res, par_res, seq_error, par_error);
	}
	printf("\n");

	// plot 1.2 - Linear
	int nListLinear[20] = {50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 550000, 600000, 650000, 700000, 750000, 800000, 850000, 900000, 950000, 1000000};
	printf("Plot 1.2 - Linear\n");
	printf("n       pi_sequential   pi_parallel     error_seq       error_par\n");
	for (int i = 0; i < 20; i++){
        double seq_res = sequentialTest(nListLinear[i]);
		double seq_error = fabs(seq_res - PI) / PI;
        double par_res = parallelTestReduce(nListLinear[i], P_DEFAULT);
		double par_error = fabs(par_res - PI) / PI;
		printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", nListLinear[i], seq_res, par_res, seq_error, par_error);
	}
	printf("\n");


    // int n = 10000000;
    // double start = omp_get_wtime();
    // parallelTestCritical(n, 1);
    // double elapsed = omp_get_wtime() - start;
    // double seq_time = elapsed;
    double start, elapsed, seq_time;

	// plot 2.1
	printf("Plot 2.1 - Critical\n");
	printf("p       pi              error           time (s)        speedup         efficiency (%%)\n");
	for (int p = 1; p <= pmax; p++){
        start = omp_get_wtime();
		double pi = parallelTestCritical(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}
    printf("\n");

    printf("Plot 2.2 - Atomic\n");
	printf("p       pi              error           time (s)        speedup         efficiency (%%)\n");
	for (int p = 1; p <= pmax; p++){
        start = omp_get_wtime();
		double pi = parallelTestAtomic(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}
    printf("\n");

    printf("Plot 2.3 - Reduction\n");
	printf("p       pi              error           time (s)        speedup         efficiency (%%)\n");
	for (int p = 1; p <= pmax; p++){
        start = omp_get_wtime();
		double pi = parallelTestReduce(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}
    printf("\n");

	return;
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int pmax = atoi(argv[2]);
    // printf("critical: %.15lf\n", parallelTestCritical(1000000, 10));
    // printf("atomic: %.15lf\n", parallelTestAtomic(1000000, 10));
    // printf("reduce: %.15lf\n", parallelTestReduce(1000000, 10));
    driver(n, pmax);
    return 0;
}