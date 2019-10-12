#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include "borg.h"

#define PI 3.14159265358979323846

int nvars = 4;
int nobjs = 3;

//Compiles with (e.g.) gcc dtlz2_single_processor.c borg.c mt19937ar.c -lm -o dtlz2_sp.exe

void dtlz2(double* vars, double* objs, double* consts) {
	int i;
	int j;
	int k = nvars - nobjs + 1;
	double g = 0.0;

	for (i=nvars-k; i<nvars; i++) {
		g += pow(vars[i] - 0.5, 2.0);
	}

	for (i=0; i<nobjs; i++) {
		objs[i] = 1.0 + g;

		for (j=0; j<nobjs-i-1; j++) {
			objs[i] *= cos(0.5*PI*vars[j]);
		}

		if (i != 0) {
			objs[i] *= sin(0.5*PI*vars[nobjs-i-1]);
		}
	}
}
// Compiles with : gcc -o dtlz2_sp.exe dtlz2_single_processor.c borg.c mt19937ar.c -lm
int main(int argc, char* argv[]) {
	int i, j;
	int rank;
	char runtime[256];
    int NFE = 10000;
	// Define the problem.  Problems are defined the same way as the
	// serial example (see dtlz2_serial.c).
	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, 0, dtlz2);

	for (j=0; j<nvars; j++) {
		BORG_Problem_set_bounds(problem, j, 0.0, 1.0);
	}

	for (j=0; j<nobjs; j++) {
		BORG_Problem_set_epsilon(problem, j, 0.01);
	}
    char buffer[100];
	// When running experiments, we want to run the algorithm multiple
	// times and average the results.
	for (i=0; i<20; i++) {


		// Seed the random number generator.
		BORG_Random_seed(37*i*(rank+1));

		// Run the Borg MOEA on the problem.
		BORG_Archive result = BORG_Algorithm_run(problem,NFE);

		// Print the Pareto optimal solutions to file.
		if (result != NULL) {
            sprintf(buffer,"/mnt/c/Users/Barney/Documents/borg/hyper_regret/%d.out",i+1);
            FILE *fid = fopen(buffer,"w");
			BORG_Archive_print(result, fid);
			BORG_Archive_destroy(result);
            fclose(fid);
		}
	}

	// Shutdown the parallel processes and exit.
	BORG_Problem_destroy(problem);
	return EXIT_SUCCESS;
}

