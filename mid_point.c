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

double originFunc(double x) {
    return 4.0 / tan(x);
}

double truth(void) {
    return originFunc(X_MAX) - originFunc(X_MIN);
}

// double sequentialTest(int n) {
// 	int sum = 0;
// 	for (int i = 0; i < n; i++) {
// 		sum += insideQuaterCircle(getRandom(), getRandom());
// 	}
// 	double pi = sum / (double) n * 4.0;
// 	// printf("[Sequential] The estimated pi = %lf\n", pi);
// 	return pi;
// }

int main(int argc, char **argv) {
    printf("%d %s\n", argc, argv[0]);
    printf("%lf\n", truth());
    return 0;
}