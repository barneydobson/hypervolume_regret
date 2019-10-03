#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <dirent.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "fex.h"

typedef struct weightEntry *node;
//Here we re-evaluate the weights on revalT (as opposed to that used for optimization) - this also makes these 'validation with framing used for optimization, rather than calibration during optimization'

struct weightEntry{
	double *weights;
	double *objectives;
	struct weightEntry *next;
	struct weightEntry *head;
};


node createEntry(){
	node temp;
	temp = (node)malloc(sizeof(struct weightEntry));
	temp->next = NULL;
	return temp;
}

node addNode(node head, double *weights, double *objectives){
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
		while(p->next != NULL){
			p = p->next;
		}
		p->next = temp;
	}
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
void printWeights(node head, int NW){
node p = head;
int i;
while (p != NULL)
{
	for (i=0;i<NW;i++)
	{
		printf("%f ", p->weights[i]);
	}
	p = p->next;
	printf("\n");
}

}
node runtimeToFinalWeights(char *fid,int numWeights, int numObjectives, node myWeights)
{
	
	int i;
	FILE *fh = fopen(fid,"r");
	char *strbuf = calloc(5000,sizeof(char));
	int offset;
	int total_offset =0;
	double temp;
	while (fgets(strbuf,4999,fh) != NULL)
	{
		if (strbuf[0] != '/' && strbuf[0] != '#')
		{
			double *weights = calloc(numWeights,sizeof(double));
								
			for (i=0;i<numWeights;i++)
			{
				sscanf(strbuf,"%lf%n",&temp,&offset);
				weights[i] = temp;
				strbuf += offset;
				total_offset += offset;
			}
			double *objectives = calloc(numObjectives,sizeof(double));
			for (i = 0;i< numObjectives ;i++)
			{
				sscanf(strbuf,"%lf%n",&temp,&offset);
				objectives[i] = temp;
				strbuf += offset;
				total_offset += offset;
			}
			strbuf = strbuf - total_offset;
			total_offset = 0;
			myWeights = addNode(myWeights,weights,objectives);
		}
		else if (strbuf[0] == '/')
		{
			myWeights = resetNode(myWeights);
		}
	}
	return myWeights;
}


int main(int argc, char* argv[])
{
	int i, j, k,l, size, rank, currentFrame = -1,currentSimType = -1;
	int masterSeed = atoi(argv[1]);
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int numFramingsForEval = 3, numFramings = 4;
	int coopIndex[4] = {1,1,1,1};
	int evalIndex[4] = {1,1,1,0};
	int framingSpecifier[4] = {1,2,3,4};
	struct fexEnsemble myEnsemble;
	int numCoopObjs = 4, numC2Objs = 2, numC1Objs = 2; //Should really do better than this - potentially move the objective options to netoptions? a lot of work so I'll leave this for now
	int objIndex[4] = {0,2,1,3}; //Same as above - not a very elegant solution (note this index describes how placing c2 objs next to c1 objs relates to the order of coop objs)
	int revalT = 30000; 
	int enSize = 400;
	struct fexEnsemble dummyEnsemble;
	dummyEnsemble.ensembleSize = 1;
	dummyEnsemble.ensembleT = 1;
	dummyEnsemble.baseSeed = 1;
	dummyEnsemble.percentile = 1;
	dummyEnsemble.simulationType = 0;
	createEnsemble(&dummyEnsemble,1);
	
	char *fid = calloc(256,sizeof(char));
	sprintf(fid,"reval_seed%d_ensembleSize%d_T%d.reeval",masterSeed,enSize,revalT);; //Consider making this .runtime - also consider a more informative filename	
	DIR *dir;
	struct dirent *ent;
	char *fStr = calloc(256,sizeof(char));
	sprintf(fStr,"optimization_seed%d_framing",1);
	
	node *allWeights = calloc(numFramings,sizeof(node));
	for (k = 0; k < numFramings;k++)
	{
		node weights = NULL;
		node c2Weights = NULL;
		node c1Weights = NULL;
		if ((dir = opendir (".")) != NULL)
		{
			while ((ent = readdir (dir)) != NULL) 
			{
				char *tempstr = malloc((strlen(fStr)+1)*sizeof(char));
				for (i = 0; i < strlen(fStr);i++)
				{
					tempstr[i] = ent->d_name[i];
				}
				tempstr[i] = '\0';
				if (strcmp(fStr,tempstr) == 0)
				{
					currentFrame = ent->d_name[strlen(fStr)] - '0'; //note, this will only work with fewer than 9 or 10 framings
				}
				else
				{
				 currentFrame=-1;
				 }
				free(tempstr);
				if (currentFrame == framingSpecifier[k])
				{
					if (coopIndex[k] == 1)
					{
						weights = runtimeToFinalWeights(ent->d_name,dummyEnsemble.ensemble[0]->cooperativeNet->numWeights,numCoopObjs,weights);
						
						if (currentFrame != 4)
						{
							struct fexEnsemble tempEnsemble;
							tempEnsemble.ensembleSize = 1;
							tempEnsemble.ensembleT = revalT;
							tempEnsemble.baseSeed = masterSeed + 40;
							tempEnsemble.percentile = 1;
							tempEnsemble.simulationType =0;
							createEnsemble(&tempEnsemble,framingSpecifier[k]);
							srand(tempEnsemble.baseSeed);
							for (j = 0; j < 1; j++)
							{
								tempEnsemble.ensemble[j]->seed = rand();
							}
							node temp = weights;
							i = 0;
							while (temp != NULL)
							{
								if (i%size == rank)
								{
									evaluateEnsembleMean(temp->weights,temp->objectives,0,&tempEnsemble);
									if (rank != 0)
									{
										MPI_Send(temp->objectives, numCoopObjs,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
									}
								}
								if (rank == 0)
								{
									if (i%size != rank)
									{
										MPI_Recv(temp->objectives, numCoopObjs,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
									}
									
								}
								temp = temp->next;
								i = i + 1;
							}	
							destroyEnsemble(&tempEnsemble);
						}		
						else
						{
							int framingsIncluded[3] = {1,2,3};
							int en2 = 1;
							struct fexEnsemble tempEnsemble = createMultiFramingEnsemble(numFramingsForEval, framingsIncluded, en2, revalT, masterSeed + 40, 50, 0);
							srand(tempEnsemble.baseSeed);
							for (j = 0; j < 1; j++)
							{
								tempEnsemble.ensemble[j]->seed = rand();
							}
							node temp = weights;
							i = 0;
							while (temp != NULL)
							{
								if (i%size == rank)
								{
									evaluateEnsembleMean(temp->weights,temp->objectives,0,&tempEnsemble);
									if (rank != 0)
									{
										MPI_Send(temp->objectives, numCoopObjs,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
									}
								}
								if (rank == 0)
								{
									if (i%size != rank)
									{
										MPI_Recv(temp->objectives, numCoopObjs,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
									}
									
								}
								temp = temp->next;
								i = i + 1;
							}	
							destroyEnsemble(&tempEnsemble);
						}
					}
					else
					{
						currentSimType = ent->d_name[strlen(fStr)+strlen("_simType")+1] - '0';
						switch (currentSimType){
							case 1:
								c2Weights = runtimeToFinalWeights(ent->d_name,dummyEnsemble.ensemble[0]->c2Net->numWeights,numC2Objs,c2Weights);
								
								struct fexEnsemble tempC2Ensemble;
								tempC2Ensemble.ensembleSize = enSize;
								tempC2Ensemble.ensembleT = revalT;
								tempC2Ensemble.baseSeed = masterSeed + 40;
								tempC2Ensemble.percentile = 1;
								tempC2Ensemble.simulationType = 1;
								createEnsemble(&tempC2Ensemble,framingSpecifier[k]);
								srand(tempC2Ensemble.baseSeed);
								for (j = 0; j < enSize; j++)
								{
									tempC2Ensemble.ensemble[j]->seed = rand();
								}
								node c2Temp = c2Weights;
								i = 0;
								while (c2Temp != NULL)
								{
									if (i%size == rank)
									{
										evaluateEnsembleMean(c2Temp->weights,c2Temp->objectives,0,&tempC2Ensemble);
										if (rank != 0)
										{
											MPI_Send(c2Temp->objectives, numC2Objs,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
										}
									}
									if (rank == 0)
									{
										if (i%size != rank)
										{
											MPI_Recv(c2Temp->objectives, numC2Objs,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
										}
										
									}
									c2Temp = c2Temp->next;
									i = i + 1;
								}	
								destroyEnsemble(&tempC2Ensemble);
								break;
							case 2:
								c1Weights = runtimeToFinalWeights(ent->d_name,dummyEnsemble.ensemble[0]->c1Net->numWeights,numC1Objs,c1Weights);
								
								struct fexEnsemble tempC1Ensemble;
								tempC1Ensemble.ensembleSize = enSize;
								tempC1Ensemble.ensembleT = revalT;
								tempC1Ensemble.baseSeed = masterSeed + 40;
								tempC1Ensemble.percentile = 1;
								tempC1Ensemble.simulationType = 2;
								createEnsemble(&tempC1Ensemble,framingSpecifier[k]);
								srand(tempC1Ensemble.baseSeed);
								for (j = 0; j < enSize; j++)
								{
									tempC1Ensemble.ensemble[j]->seed = rand();
								}
								
								node c1Temp = c1Weights;
								i = 0;
								while (c1Temp != NULL)
								{
									if (i%size == rank)
									{
										evaluateEnsembleMean(c1Temp->weights,c1Temp->objectives,0,&tempC1Ensemble);
										if (rank != 0)
										{
											MPI_Send(c1Temp->objectives, numC1Objs,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
										}
									}
									if (rank == 0)
									{
										if (i%size != rank)
										{
											MPI_Recv(c1Temp->objectives, numC1Objs,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
										}
										
									}
									c1Temp = c1Temp->next;
									i = i + 1;
								}	
								destroyEnsemble(&tempC1Ensemble);
								
								break;
						}
					}
				}				
			}
			
		}
		if (coopIndex[k] == 0){
			//Weight sets must be combined
			node tempC2 = c2Weights;
			node tempC1 = c1Weights;
			while (tempC2 != NULL)
			{
				while (tempC1 != NULL)
				{					
					double *holderWeights = calloc(dummyEnsemble.ensemble[0]->c2Net->numWeights+dummyEnsemble.ensemble[0]->c1Net->numWeights,sizeof(double));
					for (j = 0; j < dummyEnsemble.ensemble[0]->c2Net->numWeights; j++)
					{
						holderWeights[j] = tempC2->weights[j];
					}
					
					for (i = 0; i < dummyEnsemble.ensemble[0]->c1Net->numWeights;i++)
					{
						holderWeights[i+j] = tempC1->weights[i];
						
					}
					
					double *holderObjectives = calloc(numCoopObjs,sizeof(double));
					for (j=0; j < numC2Objs;j++)
					{
						holderObjectives[objIndex[j]] = tempC2->objectives[j];
					}
					for (i = 0; i< numC1Objs ; i++)
					{
						holderObjectives[objIndex[j+i]] = tempC1->objectives[i];
					}
					weights = addNode(weights,holderWeights,holderObjectives);
					tempC1 = tempC1->next;
				}
				tempC1 = c1Weights;
				tempC2 = tempC2->next;
			}
		}
		closedir(dir);
		allWeights[k] = weights;
	}
	FILE *fres;
	if (rank == 0)
	{
		fres = fopen(fid,"w");
		fprintf(fres,"//Welcome to reval, today we are using %d processors\n//each row contains %d values, half for when the objective is evaluated on itself, and half when evaluted on the framing specified\n",size,numCoopObjs);
	}
	
	struct timeval tv1, tv2;
	struct rusage ru_before, ru_after;
	

	double totaltime = 0;
	for (l=0; l < numFramings; l++)
	{
		
		if (evalIndex[l] == 1)
		{
			for (k = 0; k < numFramings; k++)
			{	
				getrusage (RUSAGE_SELF, &ru_before);
				if (rank == 0)
				{
					fprintf(fres,"//Framing %d Evaluated on Framing %d\n",framingSpecifier[k],framingSpecifier[l]);
					fflush(fres);
				}
				// fprintf(fres,"\\Taking Weights from file: %s\n", -- nice idea... could add to nodes, also calibration objectives
				myEnsemble.ensembleSize = enSize;
				myEnsemble.ensembleT = revalT;
				myEnsemble.baseSeed = masterSeed + 20;
				myEnsemble.percentile = 1;
				if (coopIndex[k] == 0)
				{
					myEnsemble.simulationType = 3;
				}
				else
				{
					myEnsemble.simulationType = 0;
				}
				createEnsemble(&myEnsemble, framingSpecifier[l]); //if this takes any time then just move out an manually change each ensemble member's simType to 0/3 (fex.c/changeFexSimType)
				srand(myEnsemble.baseSeed);
				for (j = 0; j < enSize; j++)
				{
					myEnsemble.ensemble[j]->seed = rand();
				}
				node p = allWeights[k];
		
				i = 0;
				while (p != NULL){
					double *objectives = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
					if (i%size == rank)
					{
						evaluateEnsembleMean(p->weights,objectives,0,&myEnsemble);
						if (rank != 0)
						{
							MPI_Send(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
							
						}
					}
					if (rank == 0)
					{
						if (i%size != rank)
						{
							MPI_Recv(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,i%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						}
						// printf("I completed sim, objectives[0] = %f\n", objectives[0]);
						for (j=0;j<myEnsemble.ensemble[0]->numObjectives;j++)
						{
							fprintf(fres,"%f ",p->objectives[j]);
						}
						for (j=0;j<myEnsemble.ensemble[0]->numObjectives;j++)
						{
							fprintf(fres,"%f",objectives[j]);
							if (j==(myEnsemble.ensemble[0]->numObjectives-1))
							{
								fprintf(fres,"\n");
							}
							else
							{
								fprintf(fres," ");
							}
						}
					}
					free(objectives);
					i = i + 1;
					p = p->next;
				}
				
				getrusage (RUSAGE_SELF, &ru_after);
				tv1 = ru_before.ru_utime;
				tv2 = ru_after.ru_utime;
				totaltime += tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6;
				if (rank == 0)
				{
					fprintf(fres,"#\n//Total time elapsed = %f (s)\n",totaltime);
				}
				destroyEnsemble(&myEnsemble);
				
			}
		}
	}
		
	if (rank == 0)
	{
		fclose(fres);
	}
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}