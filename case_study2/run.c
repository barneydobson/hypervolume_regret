#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "borgms.h"
#include "fex.h"
typedef struct weightEntry *node;

struct weightEntry{
	double *weights;
	double *objectives;
	struct weightEntry *next;
	struct weightEntry *previous;
	struct weightEntry *head;
};

node removeEntry(node p)
{
	node next = NULL;
	if (p != NULL)
	{
		if (p->next != NULL)
		{
			if (p->previous != NULL)
			{
				p->next->previous = p->previous;
				p->previous->next = p->next;
				next = p->next;
			}
			else
			{
				p->next->previous = NULL;
				node p_ = p->next;
				next = p->next;
				while (p_ != NULL)
				{
					p_->head = p->next;
					p_ = p_->next;
				}
			}
		}
		else
		{
			p->previous->next = NULL;
		}
		
	}
	
	free(p->weights);
	free(p->objectives);
	free(p);
	return next;
}
node createEntry(){
	node temp;
	temp = (node)malloc(sizeof(struct weightEntry));
	temp->next = NULL;
	temp->previous = NULL;
	temp->head = NULL;
	return temp;
}

node addEntry(node head, double *weights, double *objectives){
	node temp, p;
	temp = createEntry();
	temp->weights = weights;
	temp->objectives = objectives;
	if (head == NULL){
		head = temp;
	}
	else
	{
		p = head;
		while (p->next != NULL){
			p = p->next;
		}
		p->next = temp;
		temp->previous = p;
	}
	temp->head = head;
	return head;
}

node resetNode(node head){
	node current = head;
	node next;
	while (current != NULL) {
		next = current->next;
		free(current->weights);
		free(current->objectives);
		free(current);
		current = next;
	}
	return current;
}
int main(int argc, char* argv[]) {
	int masterSeed = atoi(argv[1]);
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize); //worldSize and worldRank are globals declared by borgms.c or h .. I'm hijacking them so I can use them outside of Borg as well as inside
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

	int i, j, k;
	int numWeights, numFramings = 3;
	int framingsIncluded[3] = {1,2,3};
	char framStr[100];
	// sprintf(framStr,"framings_%d_%d_%d_%d",framingsIncluded[0],framingsIncluded[1],framingsIncluded[2],framingsIncluded[3]);
	sprintf(framStr,"framing%d",numFramings+1);
	int simType = 0;	
	//-------------------------------------------------------------------//
	// Setup 
	//-------------------------------------------------------------------//
	struct fexEnsemble myEnsemble;
	myEnsemble.ensembleSize = 1;
	myEnsemble.ensembleT = 2;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 1;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = c1 only, 4 = un-cooperative
	createEnsemble(&myEnsemble, framingsIncluded[0]);	
	numWeights = myEnsemble.ensemble[0]->cooperativeNet->numWeights;
	int numObjectives = myEnsemble.ensemble[0]->numObjectives;
	double *epsilons = calloc(numObjectives,sizeof(double));
	double *epsilonUB = calloc(numObjectives,sizeof(double));
	double *epsilonLB = calloc(numObjectives,sizeof(double));
	char epstrIni[100];
	
	epsilons[0] = 0.11; // 
	epsilons[1] = 0.17; // 
	epsilons[2] = 1.3; 
	epsilons[3] = 9; 
	
	epsilonUB[0] = epsilons[0]*30;
	epsilonUB[1] = epsilons[1]*30;
	epsilonUB[2] = epsilons[2]*30;
	epsilonUB[3] = epsilons[3]*30;
	
	epsilonLB[0] = epsilons[0]/10;
	epsilonLB[1] = epsilons[1]/10;
	epsilonLB[2] = epsilons[2]/10;
	epsilonLB[3] = epsilons[3]/10;
	
	sprintf(epstrIni,"%f_%f_%f_%f",epsilons[0],epsilons[1],epsilons[2],epsilons[3]);
			
	int NFE, NRandom, nseeds;
	BORG_Problem problem;
	BORG_Archive result;
	struct timeval tv1, tv2;
	struct rusage ru_before, ru_after;
	int ensembleSizePerFraming = 134;
	int ensembleT = 30000;
	int baseSeed = masterSeed;
	int percentile = 1;
	//-------------------------------------------------------------------//
	
	//-------------------------------------------------------------------//
	// Run actual optimization
	//-------------------------------------------------------------------//
	struct fexEnsemble combinedEnsemble = createMultiFramingEnsemble(numFramings, framingsIncluded, ensembleSizePerFraming, ensembleT, baseSeed, percentile, simType);
	NRandom = 100000;
	
	
	srand(1);
	for (i = 0; i < combinedEnsemble.ensembleSize; i++)
	{
		combinedEnsemble.ensemble[i]->seed = rand();
	}
	NFE = 100000;
	
	problem = BORG_Problem_create(numWeights, numObjectives, 0, evaluateEnsembleMean, &combinedEnsemble);
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -3, 3);
	}
	for (i = 0; i < numObjectives; i++)
	{
		BORG_Problem_set_epsilon(problem, i, epsilons[i]);
	}
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);
	char runtimeFinal[512];
	sprintf(runtimeFinal,"optimization_seed%d_%s_simType%d_ensembleSizePerFraming%d_T%d_NFE%d_epsilons%s.runtime",masterSeed,framStr,myEnsemble.simulationType,ensembleSizePerFraming,ensembleT,NFE,epstrIni);
	BORG_Algorithm_output_runtime(runtimeFinal); //create runtime (i.e. evolution of performance and weights over time
	srand(masterSeed);
	BORG_Random_seed(rand());
	result = BORG_Algorithm_ms_run(problem);
	if (result != NULL)
	{
		BORG_Archive_destroy(result);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	BORG_Problem_destroy(problem);
	MPI_Finalize();
	return EXIT_SUCCESS;
	//-------------------------------------------------------------------//
	// Run fex with random weights (to show benefits of optimization) - maybe print only the paretofront please (to epsilon?)
	//-------------------------------------------------------------------//
	char randomRuntime[256];
	FILE *fidRandomRuntime;
	getrusage (RUSAGE_SELF, &ru_before);
	node randSet = NULL;
	
	int *randSeeds = malloc(NRandom * sizeof(int));
	srand(baseSeed);
	for (i = 0;i < NRandom; i++)
	{
		randSeeds[i] = rand();
	}
	
	for (i = 0;i < NRandom; i++)
	{
		double *weights = calloc(numWeights,sizeof(double));
		double *objectives = calloc(numObjectives,sizeof(double));
		srand(randSeeds[i]);
		for (j = 0; j < numWeights; j++)
		{
			weights[j] = (6*(double)rand()/(double)RAND_MAX)-3; //random weight between -1 and 1
		}
		if (i%worldSize == worldRank) 
		{
			evaluateEnsembleMean(weights,objectives,0,&combinedEnsemble);
			if (worldRank != 0)
			{
				MPI_Send(objectives,numObjectives,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
			}
		}
		if (worldRank == 0)
		{
			if (i%worldSize != worldRank)
			{
				MPI_Recv(objectives,numObjectives,MPI_DOUBLE,i%worldSize,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			randSet = addEntry(randSet,weights,objectives);
		}
		if (worldRank != 0)
		{
			free(weights);
			free(objectives);
		}
	}
	
	getrusage (RUSAGE_SELF, &ru_after);
	tv1 = ru_before.ru_utime;
	tv2 = ru_after.ru_utime;
	if (worldRank == 0)
	{
		node p1 = randSet;
		while (p1 != NULL)
		{
			node p2 = p1->head;
			while (p2 != NULL)
			{
				int index = 0;
				for (i = 0; i < numObjectives; i++)
				{
					if (p1->objectives[i] <= p2->objectives[i])
					{
						index++;
					}
				}
				if (index == numObjectives && p1 != p2)
				{
					p2 = removeEntry(p2);
				}
				else
				{
					p2 = p2->next;
				}
			}
			randSet = p1->head;
			p1 = p1->next;
		}
		sprintf(randomRuntime,"random_weights_seed%d_%s_simType%d_ensembleSizePerFraming%d_T%d_NRandom%d.runtime",masterSeed,framStr,myEnsemble.simulationType,ensembleSizePerFraming,ensembleT,NRandom);
		fidRandomRuntime = fopen(randomRuntime,"w");
		fprintf(fidRandomRuntime,"//Have a nice day, today we are using %d processors\n",worldSize);
		fprintf(fidRandomRuntime,"// total time elapsed = %f (s)\n",tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
		fflush(fidRandomRuntime);
		p1 = randSet;
		while (randSet != NULL)
		{
			for (i = 0; i < numWeights ; i++)
			{
				fprintf(fidRandomRuntime, "%f ", randSet->weights[i]);
			}
			for (i = 0; i < numObjectives ; i++)
			{
				fprintf(fidRandomRuntime, "%f", randSet->objectives[i]);
				if (i != numObjectives - 1)
				{
					fprintf(fidRandomRuntime, " ");
				}
				else
				{
					fprintf(fidRandomRuntime, "\n");
				}
			}
			randSet = randSet->next;
		}
		fprintf(fidRandomRuntime, "#");
		fclose(fidRandomRuntime);
		p1 = resetNode(p1);
		free(p1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	destroyEnsemble(&combinedEnsemble);
	
	
	MPI_Finalize();
	return EXIT_SUCCESS;
	//-------------------------------------------------------------------//
	// Run multiseed optimization
	//-------------------------------------------------------------------//
	struct fexEnsemble multiSeedEnsemble = createMultiFramingEnsemble(numFramings, framingsIncluded, ensembleSizePerFraming, ensembleT, i+1, percentile, simType);
	NFE = 100000;
	problem = BORG_Problem_create(numWeights, numObjectives, 0, evaluateEnsembleVar, &multiSeedEnsemble);
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -3, 3);
	}

	for (i = 0; i < numObjectives; i++)
	{
		BORG_Problem_set_epsilon(problem, i, epsilons[i]);
	}
	
	
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);
	
	nseeds = 10;
	
	for (i = 0;i < nseeds;i++)
	{
		srand(baseSeed + i*worldRank);
		for (j=0; j<multiSeedEnsemble.ensembleSize;j++)
		{
			multiSeedEnsemble.ensemble[j]->seed = rand();
		}
		char runtime[512];
		sprintf(runtime,"borg_multiseed_%s_simType%d_ensembleSizePerFraming%d_T%d_NFE%d_epsilons%s_seed%d.runtime",framStr,myEnsemble.simulationType,ensembleSizePerFraming,ensembleT,NFE,epstrIni,i);
		BORG_Algorithm_output_runtime(runtime); //create runtime (i.e. evolution of performance and weights over time
		BORG_Random_seed(worldSize*i*(worldRank+1));
		result = BORG_Algorithm_ms_run(problem);
		if (result != NULL)
		{
			BORG_Archive_destroy(result);	
		}
		destroyEnsemble(&multiSeedEnsemble);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	BORG_Problem_destroy(problem);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
	//-------------------------------------------------------------------//
	
	
	
	MPI_Finalize();
	return EXIT_SUCCESS;	
	
	
	
	
	
	
	
	
	//-------------------------------------------------------------------//
	// Check multiple epsilons
	//-------------------------------------------------------------------//
	int nepsilon = 1;
	struct fexEnsemble multiEpsEnsemble = createMultiFramingEnsemble(numFramings, framingsIncluded, ensembleSizePerFraming, ensembleT, baseSeed, percentile, simType);	
	NFE = 1000;

	problem = BORG_Problem_create(numWeights, numObjectives, 0, evaluateEnsembleVar, &multiEpsEnsemble);
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -1, 1);
	}
	for (i=myEnsemble.ensemble[0]->cooperativeNet->betaInd; i< myEnsemble.ensemble[0]->cooperativeNet->betaIndMax ; i++){
		BORG_Problem_set_bounds(problem, i, 0, 1);
	}
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);
		
	for (i = 0;i < nepsilon;i++)
	{
		char epstr[100];
		sprintf(epstr,"epsilons");
		double *epsilonsRand =calloc(numObjectives,sizeof(double));
		srand(i);
		for (j = 0; j < numObjectives; j++)
		{
			epsilonsRand[j] = exp(log(epsilonLB[j])+log(epsilonUB[j]/epsilonLB[j])*(double)rand()/(double)RAND_MAX);
			sprintf(epstr + strlen(epstr),"_%0.2f",epsilonsRand[j]);
			BORG_Problem_set_epsilon(problem, j, epsilonsRand[j]);
		}
		char runtime[512];
		//no need to add seed here - no seed is used - will also need to change MATLAB plotting code
		sprintf(runtime,"epislon_test_%s_simType%d_ensembleSizePerFraming%d_T%d_NFE%d_%s.runtime",framStr,myEnsemble.simulationType,ensembleSizePerFraming,ensembleT,NFE,epstr);
		BORG_Algorithm_output_runtime(runtime); //create runtime (i.e. evolution of performance and weights over time
		result = BORG_Algorithm_ms_run(problem);
		if (result != NULL)
		{
			BORG_Archive_destroy(result);	
		}
		MPI_Barrier(MPI_COMM_WORLD);
		free(epsilonsRand);
	}
	// Free any allocated memory
	
	BORG_Problem_destroy(problem);
	destroyEnsemble(&multiEpsEnsemble);
	
//-------------------------------------------------------------------//
	// Check convergence of objectives 
	// NOTE THIS IS CURRENTLY ONLY SET UP FOR AN ENSEMBLE SIZE OF 1 - EACH TIME A DIFFERENT SET OF RANDOM WEIGHTS IS EVALUATED, THE POINTER OF THE ENSEMBLE UNDER EVALUATION IS SWITCHED TO THE NRANDOM INDEX (J)
	// IF YOU WANT TO USE A LARGER ENSEMBLESIZE THAN 1 YOU HAVE TWO OPTIONS:
	// 1) CREATE AN ARRAY OF ENSEMBLES - THIS WILL BE QUICK BUT MEMORY INTENSIVE
	// 2) CREATE AND DESTROY A NEW ENSEMBLE ON EACH ITERATION INSIDE THE J LOOP - THIS WILL BE SLOWER BUT MORE MEMORY EFFICIENT
	//-------------------------------------------------------------------//
	NRandom = worldSize;
	int T0 = 1000, TMax = 5000, Nsteps = 30;

	char convergence[256];
	FILE *fidConvergence;
	
	if (worldRank==0)
	{
		sprintf(convergence,"convergence_%s_simType%d_ensembleSizePerFraming%d_TMax%d_T0%d_Nsteps%d_NRandom%d.runtime",framStr,myEnsemble.simulationType,ensembleSizePerFraming,TMax,T0,Nsteps,NRandom);
		fidConvergence = fopen(convergence,"w");
		fprintf(fidConvergence,"//Note, this isn't your average .runtime file!\n");
		fprintf(fidConvergence,"//Each new line is the weights and objectives associated with a simulation of the same weights over simulation length increasing with logspace(T0,TMax,Nsteps) file\n");
		fprintf(fidConvergence,"//The # separates each iteration of the experiment with different weights, today we are using %d processors\n",worldSize);
	}
	
	
	double totaltime = 0;
	for (i = 0; i < Nsteps ;i++)
	{
		getrusage (RUSAGE_SELF, &ru_before);
		
		double *objectives = calloc(numObjectives,sizeof(double));
		
		for (j = 0; j < NRandom; j++)
		{
			
			
			ensembleT = pow(10,log10(T0) + (log10(TMax) - log10(T0))*i/(Nsteps-1)); //logspace
			baseSeed = (Nsteps*j+i)*ensembleSizePerFraming + 1;
			srand(j+1);
			double *weights = calloc(numWeights,sizeof(double));
			for (k = 0; k < numWeights; k++)
			{
				//You can move this into if (j%worldSize..) statement if you are not certain that random numbers generated on different processors with the same seed will be identical (tests show this is the case for blue crystal)
				weights[k] = (2*(double)rand()/(double)RAND_MAX)-1; //random weight between -1 and 1
				if (k >= myEnsemble.ensemble[0]->cooperativeNet->betaInd)
				{
					if (k< myEnsemble.ensemble[0]->cooperativeNet->betaIndMax)
					{
						weights[k] = (double)rand()/(double)RAND_MAX; //random weight between 0 and 1
					}
				}
			}
			// myEnsemble.ensemble[0]->fexT = T0 + (TMax - T0)*i/(Nsteps-1); //linspace
			
			if (j%worldSize == worldRank)
			{
				struct fexEnsemble convergeEnsemble = createMultiFramingEnsemble(numFramings, framingsIncluded, ensembleSizePerFraming, ensembleT, baseSeed, percentile, simType);
				evaluateEnsembleVar(weights,objectives,0,&convergeEnsemble);
				if (worldRank != 0)
				{
					MPI_Send(objectives,numObjectives,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
				}
				destroyEnsemble(&convergeEnsemble);
			}
			if (worldRank == 0)
			{
				if (j%worldSize != 0)
				{
					MPI_Recv(objectives,numObjectives,MPI_DOUBLE,j%worldSize,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				for (k = 0; k < numWeights; k++)
				{
					fprintf(fidConvergence,"%f ", weights[k]);
				}
				for (k = 0; k < numObjectives; k++)
				{
					fprintf(fidConvergence,"%f", objectives[k]);
					if (k == (numObjectives - 1))
					{
						fprintf(fidConvergence,"\n");
					}
					else
					{
						fprintf(fidConvergence," ");
					}
					fflush(fidConvergence);
				}
			}
						
			free(weights);
		}
		if (worldRank == 0)
		{
			getrusage (RUSAGE_SELF, &ru_after);
			tv1 = ru_before.ru_utime;
			tv2 = ru_after.ru_utime;
			totaltime += tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6;
			fprintf(fidConvergence,"# //time elapsed = %f (s)\n",totaltime);
		}
		free(objectives);
		
	}
		
	if (worldRank == 0)
	{
		fclose(fidConvergence);
	}	
	MPI_Barrier(MPI_COMM_WORLD);
	destroyEnsemble(&myEnsemble);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
		
}