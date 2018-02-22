#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

unsigned int fieldVector[4];
unsigned int fieldVectorTemp[4];
int sizeX = 10;
int sizeY = 10;
int INT_SIZE = 32;

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

int getField(int index) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector >> vectorIndex;
	return vector & 1;
}

void initField(int sizeX, int sizeY) {
	for (int i = 0; i < sizeX * sizeY; i++) {
		unsetField(i, fieldVector);
	}
}

void printField() {
	for (int i = 0; i < sizeX * sizeY; i++) {
		if (i % sizeX == 0)
			printf("\n");
		printf("%d ", getField(i));
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

int countNeighbours(int index) {
	int count = 0;

// left and right
	if (checkIndex(index, index - 1) && getField(index - 1))
		count++;
	if (checkIndex(index, index + 1) && getField(index + 1))
		count++;

// above
	int row = index - sizeX;
	if (checkIndex(row, row - 1) && getField(row - 1))
		count++;
	if (checkIndex(row, row) && getField(row))
		count++;
	if (checkIndex(row, row + 1) && getField(row + 1))
		count++;

// below
	row = index + sizeX;
	if (checkIndex(row, row - 1) && getField(row - 1))
		count++;
	if (checkIndex(row, row) && getField(row))
		count++;
	if (checkIndex(row, row + 1) && getField(row + 1))
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

void copyArray(unsigned int dest[], unsigned int src[]) {
	for (int i = 0; i < 4; i++) {
		dest[i] = src[i];
	}
}

void cycle() {
	for (int i = 0; i < sizeX * sizeY; i++) {
		int neighbours = countNeighbours(i);
		if (computeFromNeighbourCount(neighbours, getField(i))) {
			setField(i, fieldVectorTemp);
		} else {
			unsetField(i, fieldVectorTemp);
		}
	}
	copyArray(fieldVector, fieldVectorTemp);
}

int main(void) {
	initField(sizeX, sizeY);
	copyArray(fieldVectorTemp, fieldVector);

	setField(1, fieldVector);
	setField(12, fieldVector);
	setField(20, fieldVector);
	setField(21, fieldVector);
	setField(22, fieldVector);

	printField();

	struct timespec start, stop;
	for (int i = 0; i < 5; i++) {
		clock_gettime(CLOCK_MONOTONIC, &start);
		cycle();
		clock_gettime(CLOCK_MONOTONIC, &stop);
		double passed = stop.tv_nsec - start.tv_nsec;
		printf("Cycle time: %f", passed);
		printField();
	}
	return EXIT_SUCCESS;
}
