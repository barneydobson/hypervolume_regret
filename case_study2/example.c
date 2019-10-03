//compiles with: gcc -o example.exe example.c fex.c -lm
#include "fex.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
int main(){
	struct fexEnsemble myEnsemble;
	myEnsemble.ensembleSize = 10;
	myEnsemble.ensembleT = 1000;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 99;
	myEnsemble.simulationType = 0;
    createEnsemble(&myEnsemble,1);
	int i;
	double *weights = calloc(myEnsemble.ensemble[0]->cooperativeNet->numWeights,sizeof(double));
	for (i = 0; i < myEnsemble.ensemble[0]->cooperativeNet->numWeights;i++)
	{
		weights[i] = 0;
	}
	double *objective = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
	myEnsemble.ensemble[0]->printDynamics=1;
	myEnsemble.ensemble[0]->dynamicFid=calloc(256,sizeof(char));
	myEnsemble.ensemble[0]->dynamicFid="/mnt/c/Users/Barney/Documents/GitHub/hypervolume_regret/test.out";
    evaluateEnsembleMean(weights, objective, 0, &myEnsemble);
    for (i = 0; i < myEnsemble.ensemble[0]->numObjectives;i++)
	{
		printf("%f ",objective[i]);
	}
	printf("\n");
	destroyEnsemble(&myEnsemble);
	free(objective);
	free(weights);
}
