#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>

static const int sizeX = 20;
static const int sizeY = 20;
static const int CHUNKS_X = 4;
static const int RUNS_PER_THREAD = 100;

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

static inline int INDEX(int x, int y, int w) {
	return (y + 1) * (w + 2) + x + 1;
}

static int evolve(bool *original, bool *next, int w, int h) {
	int change_counter = 0;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			int i = INDEX(x, y, w);

			// count neighbors
			int neighborCount = 0;
			for (int nx = x - 1; nx <= x + 1; nx++) {
				for (int ny = y - 1; ny <= y + 1; ny++) {
					if (original[INDEX(nx, ny, w)])
						neighborCount++;
				}
			}
			if (original[i])
				neighborCount--;

			if (!original[i] && neighborCount == 3) {
				// dead cell is resurrected if it has exectly three neighbors
				next[i] = 1;
			} else if (original[i]
					&& (neighborCount == 2 || neighborCount == 3)) {
				// living cell survives on 2 or three neighbors
				next[i] = 1;
			} else {
				// everything else is dead
				next[i] = 0;
			}

			change_counter += (next[i] != original[i]);
		}
	}

	return change_counter;
}

static inline void swap_vector(bool **a, bool **b) {
	bool *t = *a;
	*a = *b;
	*b = t;
}

int main(int argc, char **argv) {
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

	int w = xend - xstart;
	int h = yend - ystart;
	MPI_Datatype borderRightRec;
	MPI_Type_create_subarray(2, (int[] ) { h + 2, w + 2 },
			(int[] ) { h + 2, 1 }, (int[] ) { 0, w + 1 }, MPI_ORDER_C,
							MPI_CHAR, &borderRightRec);
	MPI_Type_commit(&borderRightRec);

	MPI_Datatype borderLeftRec;
	MPI_Type_create_subarray(2, (int[] ) { h + 2, w + 2 },
			(int[] ) { h + 2, 1 }, (int[] ) { 0, 0 }, MPI_ORDER_C,
							MPI_CHAR, &borderLeftRec);
	MPI_Type_commit(&borderLeftRec);

	MPI_Datatype borderRightSend;
	MPI_Type_create_subarray(2, (int[] ) { h + 2, w + 2 },
			(int[] ) { h + 2, 1 }, (int[] ) { 0, w }, MPI_ORDER_C,
							MPI_CHAR, &borderRightSend);
	MPI_Type_commit(&borderRightSend);

	MPI_Datatype borderLeftSend;
	MPI_Type_create_subarray(2, (int[] ) { h + 2, w + 2 },
			(int[] ) { h + 2, 1 }, (int[] ) { 0, 1 }, MPI_ORDER_C,
							MPI_CHAR, &borderLeftSend);
	MPI_Type_commit(&borderLeftSend);

	bool *field = calloc((xend - xstart + 2) * (yend - ystart + 2),
			sizeof(bool));
	bool *prev = calloc((xend - xstart + 2) * (yend - ystart + 2),
			sizeof(bool));

	if (rank == 0) {
		prev[1 * (w + 2) + 2] = true;
		prev[2 * (w + 2) + 3] = true;
		prev[3 * (w + 2) + 1] = true;
		prev[3 * (w + 2) + 2] = true;
		prev[3 * (w + 2) + 3] = true;
	}

	writeVTK(0, rank, prev, "gol", xstart, xend, ystart, yend);

	for (int cycle = 1; cycle < RUNS_PER_THREAD; cycle++) {

		// calculate
		evolve(prev, field, w, h);

		// output
		writeVTK(cycle, rank, field, "gol", xstart, xend, ystart, yend);

		// exchange
		// top -> bottom
		for (int i = 1; i < w + 1; i++) {
			field[(h + 1) * (w + 2) + i] = field[w + 2 + i];
		}

		// bottom -> top
		for (int i = 1; i < w + 1; i++) {
			field[i] = field[h * (w + 2) + i];
		}

		MPI_Request requests[4];
		MPI_Isend(field, 1, borderRightSend, rank_right, 42, communicator,
				&requests[0]);
		MPI_Isend(field, 1, borderLeftSend, rank_left, 42, communicator,
				&requests[1]);
		MPI_Irecv(field, 1, borderRightRec, rank_right, 42, communicator,
				&requests[2]);
		MPI_Irecv(field, 1, borderLeftRec, rank_left, 42, communicator,
				&requests[3]);
		MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

		swap_vector(&field, &prev);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
