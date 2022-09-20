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

double parallelTestCritical(int n, int n_thread) {
	int i;
    double sum = 0;
    double tmp;
    #pragma omp parallel num_threads(n_thread) private(tmp)
    {
        // omp_set_num_threads(n_thread);
        int tid = omp_get_thread_num();
        double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            // #pragma omp critical
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
        // omp_set_num_threads(n_thread);
        int tid = omp_get_thread_num();
        double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            // #pragma omp critical
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
        // omp_set_num_threads(n_thread);
        int tid = omp_get_thread_num();
        double part_sum = 0;
        #pragma omp for
        for (i = 0; i < n; i++) {
            tmp = 1.0 / n * func((1.0 / n) * (i + i + 1.0) / 2.0);
            // printf("t%d i= %d add %lf to %lf\n", tid, i, tmp, sum);
            // #pragma omp critical
            part_sum += tmp;
        }
        // #pragma omp atomic
        sum += part_sum;
    }
    return sum;
}


void driver(int n, int pmax) {
	// // plot 1.1 - Power
	// int nListPow[10] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
	// printf("Plot 1.1 - Power\n");
	// printf("      n    seq error      par error\n");
	// for (int i = 0; i < 10; i++){
	// 	double seqError = fabs(sequentialTest(nListPow[i]) - PI) / PI;
    //     omp_set_num_threads(50);
	// 	double parError = fabs(parallelTestReduce(nListPow[i]) - PI) / PI;
	// 	printf("%7d    %.9lf    %.9lf\n", nListPow[i], seqError, parError);
	// }
	// printf("\n");

	// // plot 1.2 - Linear
	// int nListLinear[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
	// printf("Plot 1.2 - Linear\n");
	// printf("      n    seq error      par error\n");
	// for (int i = 0; i < 10; i++){
	// 	double seqError = fabs(sequentialTest(nListLinear[i]) - PI) / PI;
    //     omp_set_num_threads(50);
	// 	double parError = fabs(parallelTestReduce(nListLinear[i]) - PI) / PI;
	// 	printf("%7d    %.9lf    %.9lf\n", nListLinear[i], seqError, parError);
	// }
	// printf("\n");


    // int n = 10000000;
    // double tSeqStart = clock();
	// sequentialTest(n);
	// double tSeq = clock() - tSeqStart;
    // printf("seq time baseline = %lf s\n", tSeq / CLOCKS_PER_SEC);
    // double start = omp_get_wtime();
    // parallelTestCritical(n, 1);
    // double elapsed = omp_get_wtime() - start;
    // double seq_time = elapsed;
    double start, elapsed, seq_time;

	// plot 2.1
	printf("Plot 2.1 - Critical\n");
	printf("n_thread  efficiency\n");
	for (int p = 1; p <= pmax; p++){
       
		// double t0 = clock();
        start = omp_get_wtime();
		double pi = parallelTestCritical(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
		// double t1 = clock();
		// double efficiency = tSeq / (t1 - t0);
		// printf("%6d    %.15lf\n", i, efficiency);
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}

    printf("Plot 2.2 - Atomic\n");
	printf("n_thread  efficiency\n");
	for (int p = 1; p <= pmax; p++){
       
		// double t0 = clock();
        start = omp_get_wtime();
		double pi = parallelTestAtomic(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
		// double t1 = clock();
		// double efficiency = tSeq / (t1 - t0);
		// printf("%6d    %.15lf\n", i, efficiency);
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}

    printf("Plot 2.3 - Reduction\n");
	printf("n_thread  efficiency\n");
	for (int p = 1; p <= pmax; p++){
       
		// double t0 = clock();
        start = omp_get_wtime();
		double pi = parallelTestReduce(n, p);
        elapsed = omp_get_wtime() - start;
        if (p == 1) seq_time = elapsed;
		// double t1 = clock();
		// double efficiency = tSeq / (t1 - t0);
		// printf("%6d    %.15lf\n", i, efficiency);
        printf("%d\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", p, pi, fabs(pi - PI) / PI, elapsed, seq_time/elapsed, 100.0*seq_time/elapsed/p);
	}
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