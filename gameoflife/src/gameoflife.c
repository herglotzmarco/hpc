#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

static const int sizeX = 3000;
static const int sizeY = 3000;
static const int INT_SIZE = 32;
static const int MAX_THREADS = 2;
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
	for (int col = 0; col < sizeX; col++) {
		for (int row = 0; row < sizeY; row++) {
			printf("%d ", getField(row * sizeX + col, fieldVector));
		}
		printf("\n");
	}
	printf("\n-----------------------------------\n");
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

void cycle(bitvector *fieldVector, bitvector *nextFieldVector,
		int fieldVectorLength, domain *domains) {
#pragma omp parallel
	{
		printf("Thread %d starting subdomain\n", omp_get_thread_num());
		cycleSubdomain(domains[omp_get_thread_num()], fieldVector,
				nextFieldVector);
		printf("Thread %d finished subdomain. Waiting...\n",
				omp_get_thread_num());
	}
	printf("All threads finished and synchronized\n");
}

void domainDecomposition(domain *domains, int threadCount) {
	int nextRowStart = 0;
	int nextColStart = 0;

	for (int i = 0; i < threadCount; i++) {
		domains[i].rowStart = nextRowStart;
		domains[i].colStart = nextColStart;
		domains[i].colEnd = sizeX;
		if (i == threadCount - 1) {
			domains[i].rowEnd = sizeY;
			nextRowStart = domains[i].rowEnd;
		} else {
			domains[i].rowEnd = sizeY / threadCount;
			nextRowStart = domains[i].rowEnd;
		}
	}
}

void cycleAndMeasureTime(bitvector *fieldVector, bitvector *nextFieldVector,
		int fieldVectorLength, int threadCount) {
	omp_set_num_threads(threadCount);
	domain domains[threadCount];
	domainDecomposition(domains, threadCount);

	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	cycle(fieldVector, nextFieldVector, fieldVectorLength, domains);
	clock_gettime(CLOCK_MONOTONIC, &end);

	double elapsedSeconds = (end.tv_sec - start.tv_sec) * 1E9;
	double elapsedNanos = end.tv_nsec - start.tv_nsec;
	double totalElapsedNanos = elapsedSeconds + elapsedNanos;
	printf("Elapsed time during cycle with %d threads: %fms\n\n", threadCount,
			totalElapsedNanos / 1E6);
}

void writeVTK2(int id, bitvector *fieldVector, char prefix[1024], int w, int h) {
  char filename[2048];
  //int x,y;

  int offsetX=0;
  int offsetY=0;
  float deltax=1.0;
  //float deltay=1.0;
  long  nxy = w * h * sizeof(float);

  snprintf(filename, sizeof(filename), "%s%d%s", prefix,id, ".vti");
  FILE* fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
  fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

  for (int col = 0; col < sizeX; col++) {
  		for (int row = 0; row < sizeY; row++) {
  			float value = getField(row * sizeX + col, fieldVector);
  			fwrite((unsigned char*)&value, sizeof(float), 1, fp);
  		}
  		printf("\n");
  	}

  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void cycleAndMeasureTimeWithoutPrint(int fieldVectorLength,
		bitvector *fieldVector, bitvector *nextFieldVector) {
	for (int threads = 1; threads <= MAX_THREADS; threads++) {
		for (int i = 0; i < RUNS_PER_THREAD; i++) {
			cycleAndMeasureTime(fieldVector, nextFieldVector, fieldVectorLength,
					threads);
			swapArray(&fieldVector, &nextFieldVector);
		}
	}
}

void cycleAndMeasureTimeWithPrint(int fieldVectorLength, bitvector *fieldVector,
		bitvector *nextFieldVector) {
	setField(1, fieldVector);
	setField(12, fieldVector);
	setField(20, fieldVector);
	setField(21, fieldVector);
	setField(22, fieldVector);
	printField(fieldVector);

	for (int threads = 1; threads <= MAX_THREADS; threads++) {
		for (int i = 0; i < RUNS_PER_THREAD; i++) {
			cycleAndMeasureTime(fieldVector, nextFieldVector, fieldVectorLength,
					threads);
			swapArray(&fieldVector, &nextFieldVector);
			printField(fieldVector);
		}
	    writeVTK2(threads, fieldVector,"gol", sizeX, sizeY);
	}
}

int main(void) {
	int fieldVectorLength = (sizeX * sizeY / INT_SIZE) + 1;
	bitvector *fieldVector = calloc(fieldVectorLength, sizeof(bitvector));
	bitvector *nextFieldVector = calloc(fieldVectorLength, sizeof(bitvector));

	cycleAndMeasureTimeWithoutPrint(fieldVectorLength, fieldVector,
			nextFieldVector);

	free(fieldVector);
	free(nextFieldVector);
	return EXIT_SUCCESS;
}
