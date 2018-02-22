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

void cycle() {

}


int main(void) {
	initField(sizeX, sizeY);
	setField(1);
	setField(12);
	setField(20);
	setField(21);
	setField(22);
	printField();

	for (int i = 0; i < 10; i++) {
		cycle();
		printField();
	}
	return EXIT_SUCCESS;
}
