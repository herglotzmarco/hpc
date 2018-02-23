#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

static const int sizeX = 10;
static const int sizeY = 10;
static const int INT_SIZE = 32;

typedef unsigned int bitvector;

void setField(int index, bitvector *fieldVector) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector | (1 << vectorIndex);
	fieldVector[fieldIndex] = vector;
}

void unsetField(int index, bitvector *fieldVector) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector & ~(1 << vectorIndex);
	fieldVector[fieldIndex] = vector;
}

int getField(int index, bitvector *fieldVector) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector >> vectorIndex;
	return vector & 1;
}

void printField(bitvector *fieldVector) {
	for (int col = 0; col < sizeX; col++) {
		for (int row = 0; row < sizeY; row++) {
			printf("%d ", getField(row * sizeX + col, fieldVector));
		}
		printf("\n");
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

int countNeighbours(int col, int row, bitvector *fieldVector) {
	int count = 0;
	int index = row * sizeX + col;

// left and right
	if (checkIndex(index, index - 1) && getField(index - 1, fieldVector))
		count++;
	if (checkIndex(index, index + 1) && getField(index + 1, fieldVector))
		count++;

// above
	int offset = index - sizeX;
	if (checkIndex(offset, offset - 1) && getField(offset - 1, fieldVector))
		count++;
	if (checkIndex(offset, offset) && getField(offset, fieldVector))
		count++;
	if (checkIndex(offset, offset + 1) && getField(offset + 1, fieldVector))
		count++;

// below
	offset = index + sizeX;
	if (checkIndex(offset, offset - 1) && getField(offset - 1, fieldVector))
		count++;
	if (checkIndex(offset, offset) && getField(offset, fieldVector))
		count++;
	if (checkIndex(offset, offset + 1) && getField(offset + 1, fieldVector))
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

void swapArray(bitvector **dest, bitvector **src) {
	bitvector *temp = *dest;
	*dest = *src;
	*src = temp;
}

void cycle(bitvector *fieldVector, bitvector *fieldVectorTemp,
		int fieldVectorLength) {
#pragma omp parallel for
	for (int col = 0; col < sizeX; col++) {
		for (int row = 0; row < sizeY; row++) {
			int neighbours = countNeighbours(col, row, fieldVector);
			if (computeFromNeighbourCount(neighbours,
					getField(row * sizeX + col, fieldVector))) {
				setField(row * sizeX + col, fieldVectorTemp);
			} else {
				unsetField(row * sizeX + col, fieldVectorTemp);
			}
		}
	}
}

void cycleAndMeasureTime(bitvector *fieldVector, bitvector *fieldVectorTemp,
		int fieldVectorLength, int threadCount) {
	omp_set_num_threads(threadCount);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	cycle(fieldVector, fieldVectorTemp, fieldVectorLength);
	clock_gettime(CLOCK_MONOTONIC, &end);

	double elapsedSeconds = (end.tv_sec - start.tv_sec) * 1E9;
	double elapsedNanos = end.tv_nsec - start.tv_nsec;
	double totalElapsedNanos = elapsedSeconds + elapsedNanos;
	//		printf("elapsed seconds %f\n", elapsedSeconds);
	//		printf("elapsed nanos %f\n", elapsedNanos);
	//		printf("elapsed time total %f\n", totalElapsedNanos);
	printf("Elapsed time during cycle with %d threads: %fms\n", threadCount,
			totalElapsedNanos / 1E6);
}

int main(void) {
	int fieldVectorLength = (sizeX * sizeY / INT_SIZE) + 1;
	bitvector *fieldVector = calloc(fieldVectorLength, sizeof(bitvector));
	bitvector *fieldVectorTemp = calloc(fieldVectorLength, sizeof(bitvector));

	setField(1, fieldVector);
	setField(12, fieldVector);
	setField(20, fieldVector);
	setField(21, fieldVector);
	setField(22, fieldVector);

	printField(fieldVector);

	for (int i = 0; i < 5; i++) {
		cycleAndMeasureTime(fieldVector, fieldVectorTemp, fieldVectorLength, 1);
		swapArray(&fieldVector, &fieldVectorTemp);
		printField(fieldVector);
	}
	for (int i = 0; i < 5; i++) {
		cycleAndMeasureTime(fieldVector, fieldVectorTemp, fieldVectorLength, 2);
		swapArray(&fieldVector, &fieldVectorTemp);
		printField(fieldVector);
	}
	return EXIT_SUCCESS;
}
