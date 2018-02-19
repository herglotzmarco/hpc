#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#define TRYS 5000000

static int throw() {
	double x, y;
	x = (double) rand() / (double) RAND_MAX;
	y = (double) rand() / (double) RAND_MAX;
	if ((x * x + y * y) <= 1.0)
		return 1;

	return 0;
}

// ohne reduction
/*
 int main(int argc, char **argv) {
 int globalCount = 0, globalSamples = TRYS;

 #pragma omp parallel for
 for (int i = 0; i < globalSamples; ++i) {
 // atomic nur für addition selbst, throw() kann parallelisiert sein
 int r = throw();
 #pragma omp atomic
 globalCount += r;
 }

 double pi = 4.0 * (double) globalCount / (double) (globalSamples);

 printf("pi is %.9lf\n", pi);

 return 0;
 }
 */

// mit reduction
/*
 int main(int argc, char **argv) {
 int globalCount = 0, globalSamples = TRYS;

 #pragma omp parallel for reduction( +:globalCount)
 for (int i = 0; i < globalSamples; ++i) {
 globalCount += throw();
 }
 double pi = 4.0 * (double) globalCount / (double) (globalSamples);

 printf("pi is %.9lf\n", pi);

 return 0;
 }
 */

// mit trefferausgabe
int main(int argc, char **argv) {
	int globalCount = 0, globalSamples = TRYS;
	long int thread_count = 6;
	if (argc > 1)
		thread_count = strtol(argv[1], NULL, 0);
	if (thread_count < 1)
		thread_count = 6;
	omp_set_num_threads(thread_count);

#pragma omp parallel reduction(+:globalCount)
	{
#pragma omp for
		for (int i = 0; i < globalSamples; ++i) {
			globalCount += throw();
		}
		printf("Thread %d: Trefferanzahl: %d\n", omp_get_thread_num(),
				globalCount);
	}
	double pi = 4.0 * (double) globalCount / (double) (globalSamples);

	printf("pi is %.9lf\n", pi);

	return 0;
}
