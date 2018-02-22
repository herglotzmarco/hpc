#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

static const int sizeX = 5000;
static const int sizeY = 1500;
static const int INT_SIZE = 32;

void setField(int index, unsigned int fieldVector[]) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector | (1 << vectorIndex);
	fieldVector[fieldIndex] = vector;
}

void unsetField(int index, unsigned int fieldVector[]) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector & ~(1 << vectorIndex);
	fieldVector[fieldIndex] = vector;
}

int getField(int index, unsigned int fieldVector[]) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector >> vectorIndex;
	return vector & 1;
}

void initField(int sizeX, int sizeY, unsigned int fieldVector[]) {
	for (int i = 0; i < sizeX * sizeY; i++) {
		unsetField(i, fieldVector);
	}
}

void printField(unsigned int fieldVector[]) {
	for (int i = 0; i < sizeX * sizeY; i++) {
		if (i % sizeX == 0)
			printf("\n");
		printf("%d ", getField(i, fieldVector));
	}
	printf("\n\n-----------------------------------\n");
}

int checkIndex(int current, int i) {
	int bounds = i >= 0 && i < sizeX * sizeY;
	int rowBoundsRight = current % sizeX == sizeX - 1;
	if (rowBoundsRight) {
		rowBoundsRight = rowBoundsRight && i != current + 1;
	} else {
		rowBoundsRight = 1;
	}

	int rowBoundsLeft = current % sizeX == 0;
	if (rowBoundsLeft) {
		rowBoundsLeft = rowBoundsLeft && i != current - 1;
	} else {
		rowBoundsLeft = 1;
	}

	return bounds && rowBoundsLeft && rowBoundsRight;
}

int countNeighbours(int index, unsigned int fieldVector[]) {
	int count = 0;

// left and right
	if (checkIndex(index, index - 1) && getField(index - 1, fieldVector))
		count++;
	if (checkIndex(index, index + 1) && getField(index + 1, fieldVector))
		count++;

// above
	int row = index - sizeX;
	if (checkIndex(row, row - 1) && getField(row - 1, fieldVector))
		count++;
	if (checkIndex(row, row) && getField(row, fieldVector))
		count++;
	if (checkIndex(row, row + 1) && getField(row + 1, fieldVector))
		count++;

// below
	row = index + sizeX;
	if (checkIndex(row, row - 1) && getField(row - 1, fieldVector))
		count++;
	if (checkIndex(row, row) && getField(row, fieldVector))
		count++;
	if (checkIndex(row, row + 1) && getField(row + 1, fieldVector))
		count++;

	return count;
}

int computeFromNeighbourCount(int neighbours, int current) {
	// dead cell
	if (!current && neighbours == 3)
		return 1; // come to life
	if (!current)
		return 0; // stay dead

	// living cells
	if (neighbours < 2)
		return 0; // die of starvation
	if (neighbours >= 4)
		return 0; // die of overpopulation
	return 1; // else stay alive
}

void copyArray(unsigned int dest[], unsigned int src[], int fieldVectorLength) {
	for (int i = 0; i < fieldVectorLength; i++) {
		dest[i] = src[i];
	}
}

void cycle(unsigned int fieldVector[], unsigned int fieldVectorTemp[],
		int fieldVectorLength) {
	for (int i = 0; i < sizeX * sizeY; i++) {
		int neighbours = countNeighbours(i, fieldVector);
		if (computeFromNeighbourCount(neighbours, getField(i, fieldVector))) {
			setField(i, fieldVectorTemp);
		} else {
			unsetField(i, fieldVectorTemp);
		}
	}
	copyArray(fieldVector, fieldVectorTemp, fieldVectorLength);
}

int main(void) {
	int fieldVectorLength = (sizeX * sizeY / INT_SIZE) + 1;
	unsigned int fieldVector[fieldVectorLength];
	unsigned int fieldVectorTemp[fieldVectorLength];
	initField(sizeX, sizeY, fieldVector);
	copyArray(fieldVectorTemp, fieldVector, fieldVectorLength);

	setField(1, fieldVector);
	setField(12, fieldVector);
	setField(20, fieldVector);
	setField(21, fieldVector);
	setField(22, fieldVector);

	struct timespec start, end;
	for (int i = 0; i < 10; i++) {
		clock_gettime(CLOCK_MONOTONIC, &start);
		cycle(fieldVector, fieldVectorTemp, fieldVectorLength);
		clock_gettime(CLOCK_MONOTONIC, &end);

		double elapsedSeconds = (end.tv_sec - start.tv_sec) * 1E9;
		double elapsedNanos = end.tv_nsec - start.tv_nsec;
		double totalElapsedNanos = elapsedSeconds + elapsedNanos;

//		printf("elapsed seconds %f\n", elapsedSeconds);
//		printf("elapsed nanos %f\n", elapsedNanos);
//		printf("elapsed time total %f\n", totalElapsedNanos);

		printf("Elapsed time during cycle: %fms\n", totalElapsedNanos / 1E6);
	}
	return EXIT_SUCCESS;
}
