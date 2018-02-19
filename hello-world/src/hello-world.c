#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

int main(int argc, char **argv) {

// Aufgabe 1 b)
	/*
	 #pragma omp parallel for num_threads(4)
	 for (int i = 0; i < 4; i++) {
	 printf("Hello World from thread %d of %d\n", omp_get_thread_num(),
	 omp_get_num_threads());
	 }
	 return 0;
	 }
	 */

	/* Aufgabe 1 c)
	 * Reihenfolge ist zufällig
	 */

// Aufgabe 1 d)
	/*
	 #pragma omp parallel sections num_threads(4)
	 {
	 printf("Hallo Welt from thread %d of %d\n", omp_get_thread_num(),
	 omp_get_num_threads());
	 #pragma omp section
	 {
	 sleep(1);
	 printf("Hello World from thread %d of %d\n", omp_get_thread_num(),
	 omp_get_num_threads());
	 }
	 #pragma omp section
	 {
	 sleep(1);
	 printf("Bonjour tout le monde from thread %d of %d\n",
	 omp_get_thread_num(), omp_get_num_threads());
	 }
	 #pragma omp section
	 {
	 sleep(1);
	 printf("Hej varlden from thread %d of %d\n", omp_get_thread_num(),
	 omp_get_num_threads());
	 }
	 #pragma omp section
	 {
	 sleep(1);
	 printf("Hola mundo from thread %d of %d\n", omp_get_thread_num(),
	 omp_get_num_threads());
	 }
	 }
	 */
}
