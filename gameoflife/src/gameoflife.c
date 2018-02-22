#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

unsigned int fieldVector[4];
int sizeX = 10;
int sizeY = 10;
int INT_SIZE = 32;

void setField(int index) {
	int fieldIndex = index / INT_SIZE;
	int vectorIndex = index % INT_SIZE;
	int vector = fieldVector[fieldIndex];
	vector = vector | (1 << vectorIndex);
	fieldVector[fieldIndex] = vector;
}

void unsetField(int index) {
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
		unsetField(i);
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

void cycle() {
	for (int i = 0; i < sizeX * sizeY; i++) {
		int neighbours = countNeighbours(i);
		if (computeFromNeighbourCount(neighbours, getField(i))) {
			setField(i);
		} else {
			unsetField(i);
		}
	}
}

int main(void) {
	initField(sizeX, sizeY);
	setField(1);
	setField(12);
	setField(20);
	setField(21);
	setField(22);
	printField();

	for (int i = 0; i < 1; i++) {
		cycle();
		printField();
	}
	return EXIT_SUCCESS;
}
