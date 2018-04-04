#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>

static const int sizeX = 20;
static const int sizeY = 20;
static const int CHUNKS_X = 4;
static const int CHUNKS_Y = 1;
static const int INT_SIZE = 32;
static const int RUNS_PER_THREAD = 10;

typedef unsigned int bitvector;

struct domain {
	int rowStart;
	int rowEnd;
	int colStart;
	int colEnd;
};
typedef struct domain domain;

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
	for (int row = 0; row < sizeY; row++) {
		for (int col = 0; col < sizeX; col++) {
			printf("%d ", getField(row * sizeX + col, fieldVector));
		}
		printf("\n");
	}
	printf("\n-----------------------------------\n");
}

void writeVTK(int t, int thread, bool *field, char prefix[1024], int xleft,
		int xright, int ytop, int ybottom) {
	char name[1024] = "\0";
	sprintf(name, "%s_%02d_%04d.vtk", prefix, thread, t);
	FILE* outfile = fopen(name, "w");

	/*Write vtk header */
	fprintf(outfile, "# vtk DataFile Version 3.0\n");
	fprintf(outfile, "frame %d\n", t);
	fprintf(outfile, "BINARY\n");
	fprintf(outfile, "DATASET STRUCTURED_POINTS\n");
	fprintf(outfile, "DIMENSIONS %d %d %d \n", xright - xleft, ybottom - ytop,
			1);
	// paraview weirdness - this gets rid of spaces between the chunks
	fprintf(outfile, "SPACING %f %f 1.0\n",
			((double) xright - xleft) / ((double) xright - xleft - 1),
			((double) ybottom - ytop) / ((double) ybottom - ytop - 1));
	fprintf(outfile, "ORIGIN %d %d 0\n", xleft, ytop);
	fprintf(outfile, "POINT_DATA %d\n", (xright - xleft) * (ybottom - ytop));
	fprintf(outfile, "SCALARS data unsigned_char 1\n");
	fprintf(outfile, "LOOKUP_TABLE default\n");

//    static_assert(sizeof(bool) == sizeof(char), "bool is char");

	int w = xright - xleft;
	int h = ybottom - ytop;
	for (int row = 1; row < h + 1; row++) {
		for (int col = 1; col < w + 1; col++) {
			bool value = field[row * (w + 2) + col];
			fwrite((unsigned char*) &value, sizeof(bool), 1, outfile);
		}
	}
	fclose(outfile);
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

void cycleSubdomain(domain d, bitvector *fieldVector,
		bitvector *nextFieldVector) {
	for (int row = d.rowStart; row < d.rowEnd; row++) {
		for (int col = d.colStart; col < d.colEnd; col++) {
			int neighbours = countNeighbours(col, row, fieldVector);
			if (computeFromNeighbourCount(neighbours,
					getField(row * sizeX + col, fieldVector))) {
				setField(row * sizeX + col, nextFieldVector);
			} else {
				unsetField(row * sizeX + col, nextFieldVector);
			}
		}
	}
}



void cycle(int cycleNum, bitvector *fieldVector, bitvector *nextFieldVector,
		int fieldVectorLength, domain *domains) {
	for (int i = 0; i < CHUNKS_X * CHUNKS_Y; i++) {
		printf(
				"Domain %d: rowStart: %d, rowEnd: %d, colStart: %d, colEnd: %d\n",
				i, domains[i].rowStart, domains[i].rowEnd, domains[i].colStart,
				domains[i].colEnd);
		cycleSubdomain(domains[i], fieldVector, nextFieldVector);
//		writeVTK2(cycleNum, i, fieldVector, "gol", sizeX, sizeY, domains);
	}
	printf("All threads finished and synchronized\n");
}

void cycleAndMeasureTime(int cycleNum, bitvector *fieldVector,
		bitvector *nextFieldVector, int fieldVectorLength) {
	domain domains[CHUNKS_X * CHUNKS_Y];

//	writePVTK(cycleNum, "gol", domains);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	cycle(cycleNum, fieldVector, nextFieldVector, fieldVectorLength, domains);
	clock_gettime(CLOCK_MONOTONIC, &end);

	double elapsedSeconds = (end.tv_sec - start.tv_sec) * 1E9;
	double elapsedNanos = end.tv_nsec - start.tv_nsec;
	double totalElapsedNanos = elapsedSeconds + elapsedNanos;
	printf("Elapsed time during cycle: %fms\n\n", totalElapsedNanos / 1E6);
}

void cycleAndMeasureTimeWithoutPrint(int fieldVectorLength,
		bitvector *fieldVector, bitvector *nextFieldVector) {
	for (int i = 0; i < RUNS_PER_THREAD; i++) {
		cycleAndMeasureTime(i, fieldVector, nextFieldVector, fieldVectorLength);
		swapArray(&fieldVector, &nextFieldVector);
	}
}

void cycleAndMeasureTimeWithPrint(int fieldVectorLength, bitvector *fieldVector,
		bitvector *nextFieldVector) {
	setField(0 * sizeY + 1, fieldVector);
	setField(1 * sizeY + 2, fieldVector);
	setField(2 * sizeY + 0, fieldVector);
	setField(2 * sizeY + 1, fieldVector);
	setField(2 * sizeY + 2, fieldVector);
	printField(fieldVector);

	for (int i = 0; i < RUNS_PER_THREAD; i++) {
		cycleAndMeasureTime(i, fieldVector, nextFieldVector, fieldVectorLength);
		swapArray(&fieldVector, &nextFieldVector);
		printField(fieldVector);
	}
}

void initMPI(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int dims = { CHUNKS_X };
	int periodic = { 1 };
	MPI_Comm communicator;
	MPI_Cart_create(MPI_COMM_WORLD, 1, &dims, &periodic, 1, &communicator);

	int rank;
	MPI_Comm_rank(communicator, &rank);
	int rank_left;
	int rank_right;
	MPI_Cart_shift(communicator, 0, 1, &rank_left, &rank_right);

	int coords;
	MPI_Cart_coords(communicator, rank, 1, &coords);
	int xstart = coords * sizeX / CHUNKS_X;
	int ystart = 0;
	int xend = xstart + sizeX / CHUNKS_X;
	int yend = sizeY;

	printf(
			"My rank is %d. My neighbours are %d and %d. My coordinates are %d\n",
			rank, rank_left, rank_right, coords);
	printf(
			"I will calculate the area: xstart=%d, xend=%d, ystart=%d, yend=%d\n",
			xstart, xend, ystart, yend);

	bool *field = calloc((xend - xstart + 2) * (yend - ystart + 2),
			sizeof(bool));
	writeVTK(0, rank, field, "", xstart, xend, ystart, yend);

}

int main(int argc, char **argv) {
	int fieldVectorLength = (sizeX * sizeY / INT_SIZE) + 1;
	bitvector *fieldVector = calloc(fieldVectorLength, sizeof(bitvector));
	bitvector *nextFieldVector = calloc(fieldVectorLength, sizeof(bitvector));

	initMPI(argc, argv);
//	cycleAndMeasureTimeWithPrint(fieldVectorLength, fieldVector,
//			nextFieldVector);

	free(fieldVector);
	free(nextFieldVector);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
