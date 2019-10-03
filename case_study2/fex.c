#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
// #include "qr_solve.h"
#include "fex.h"


double doumin(double *vec,int len)
{
	//Minimum value in a vector (vec) of doubles (length = len)
	double minno = DBL_MAX;
	int i;
	for (i=0;i<len;i++)
	{
		if (minno > vec[i])
		{
			minno = vec[i];
		}
	}
	return minno;
}
double doumax(double *vec,int len)
{
	//Maximum value in a vector (vec) of doubles (length = len)
	int minno = DBL_MIN;
	int i;
	for (i=0;i<len;i++)
	{
		if (minno < vec[i])
		{
			minno = vec[i];
		}
	}
	return minno;
}
void changeFexSimType(struct fexProblem *myFex, int simulationType)
{
	//Changes a fexProblem's simulation type (0 = cooperative, 1 = c2 only, 2 = c1 only, 3 = uncooperative (weights for the two policies are joined into one weights vector with c2 policy weights first then c1 policy weights)
	myFex->simulationType = simulationType;
	if (simulationType == 0 || simulationType == 3)
	{
		myFex->numObjectives=4;//do better
		
		myFex->objectiveSpecifier[0] = 0; //C2 reliability
		myFex->objectiveSpecifier[1] = 0; //C1 reliability
		myFex->objectiveSpecifier[2] = 1; //C2 linear def
		myFex->objectiveSpecifier[3] = 0; //C2 sqr def
		myFex->objectiveSpecifier[4] = 1; //C1 linear def
		myFex->objectiveSpecifier[5] = 0; //C1 sqr def
		myFex->objectiveSpecifier[6] = 0; //netNPV
		myFex->objectiveSpecifier[7] = 0; //C2NPV
		myFex->objectiveSpecifier[8] = 0; //C1NPV
		myFex->objectiveSpecifier[9] = 0; //C1 lin gs
		myFex->objectiveSpecifier[10] = 0; //C1 sqr gs
		myFex->objectiveSpecifier[11] = 0; //net
		myFex->objectiveSpecifier[12] = 1; //C2
		myFex->objectiveSpecifier[13] = 1; //C1

	}
	else if (simulationType == 1)
	{
		myFex->numObjectives=2;
		
		myFex->objectiveSpecifier[0] = 0; 
		myFex->objectiveSpecifier[1] = 0; 
		myFex->objectiveSpecifier[2] = 1; 
		myFex->objectiveSpecifier[3] = 0; 
		myFex->objectiveSpecifier[4] = 0; 
		myFex->objectiveSpecifier[5] = 0; 
		myFex->objectiveSpecifier[6] = 0; 
		myFex->objectiveSpecifier[7] = 0; 
		myFex->objectiveSpecifier[8] = 0; 
		myFex->objectiveSpecifier[9] = 0; 
		myFex->objectiveSpecifier[10] = 0;
		myFex->objectiveSpecifier[11] = 0;
		myFex->objectiveSpecifier[12] = 1;
		myFex->objectiveSpecifier[13] = 0;
		

	}
	else if (simulationType == 2)
	{
		myFex->numObjectives=2;
		
		myFex->objectiveSpecifier[0] = 0; 
		myFex->objectiveSpecifier[1] = 0;
		myFex->objectiveSpecifier[2] = 0;
		myFex->objectiveSpecifier[3] = 0;
		myFex->objectiveSpecifier[4] = 1;
		myFex->objectiveSpecifier[5] = 0;
		myFex->objectiveSpecifier[6] = 0;
		myFex->objectiveSpecifier[7] = 0;
		myFex->objectiveSpecifier[8] = 0;
		myFex->objectiveSpecifier[9] = 0;
		myFex->objectiveSpecifier[10] = 0;
		myFex->objectiveSpecifier[11] = 0;
		myFex->objectiveSpecifier[12] = 0;
		myFex->objectiveSpecifier[13] = 1;
	
	}
	
}
int sumInt(int *vec, int T)
{
	//Sum of vector (vec) of ints (length = T)
	int i;
	int sum = 0;
	for (i=0;i<T;i++)
	{
		sum = sum + vec[i];
	}
	return sum;
}

void createFex(struct fexProblem *myFex, struct fexEnsemble *myEnsemble, int T, int seed, int simulationType)
{
	//Create a fex problem - in this function you can alter the parameters for simulation, choose what variables are used in a policy and specify the architecture of the policy network.
	//Memory is allocated here for forcing variables (although these forcing variables are populated in fex.c/populateSeededVariables)
	//fexEnsemble is only required here since it contains historical data - nothing else from it is used
	//NO NEED FOR BASESEED HERE
	int i,j;
	/// -------------- SYSTEM PARAMETERS --------------- ///
	myFex->linkCapS2D2=36.0; // daily link cap in Ml
	myFex->linkCapS1D2=45.0; // daily link cap in Ml
	myFex->linkCapR1S1=150; // daily link cap in Ml
	myFex->linkCapS1R1=600;// daily link cap in Ml
	myFex->linkCapBorehole=50;// daily link cap in Ml
	myFex->linkAnnualCapS1D2 = 14917.26;// annual link cap in Ml
	myFex->linkAnnualCapR1S1 = 13633; // annual link cap in Ml
	myFex->linkAnnualCapS2D2 = 10000;// annual link cap in Ml
	myFex->linkMonthlyCapS2D2 = 1120;// monthly link cap in Ml
	myFex->linkAnnualCapR1S1 = 12585;// annual link cap in Ml - note that this limit is ONLY on water that is later consumed by c1
	
	myFex->storCapS2=5662.21;
	myFex->storCapS1=21365;
	
	myFex->discountRate = 0.03/365;
	
	myFex->deadStorS2=0;
	myFex->deadStorS1=0;
	
	myFex->startDoyNoR1S1 = 91; //1st april -- no R1S1 between these dates
	myFex->endDoyNoR1S1 = 304; //31st october
	
	myFex->startDoyMaxR1D2 = 305; //In simType 2, this is when uS1D2 = max
	myFex->endDoyMaxR1D2 = 90; //the rest of the annual license is evenly spread outside these dates
	
	myFex->storElevationS2[0] = 1.564; //Elev = param[0]*stor^param[1]
	myFex->storElevationS2[1] = 0.3418;
	myFex->elevationAreaS2[0] = 3200;
	myFex->elevationAreaS2[1] = 1.5;
	
	myFex->storElevationS1[0] = 3.849;
	myFex->storElevationS1[1] = 0.2569; 
	myFex->elevationAreaS1[0] = 4250;
	myFex->elevationAreaS1[1] = 1.5;
	
	
	myFex->s1CompensationFixed = 1;
	myFex->s2CompensationFixed = 5;
	
	myFex->R1Thresh = 1.16*86.4;
	myFex->R2Thresh = 3.158*86.4;
		
	myFex->fexDemandD1Fixed = 200;
	
	myFex->r2r1ConversionFactor = 1/0.746; //(Flow at R2) = r2r1ConversionFactor * (Flow at R1), C2 use this conversion factor
	/// ------------------------------------------------ ///
	
	
	/// ------------- SIMULATION PARAMETERS ------------ ///
	myFex->fexInitialStorS2=5000; //Note - this may be changed in fex.c/populateSeededVariables
	myFex->fexInitialStorS1=18000; //Note - this may be changed in fex.c/populateSeededVariables
	myFex->simulationType = simulationType; //don't change this here
	myFex->printDynamics = 0; //don't change this here
	// myFex->includetimeObjectives = 0;
	myFex->fexT = T;
	//edit: LoRam 4/Apr/18
	myFex->warmup = 1000;
	myFex->seed = seed;
	//endedit

	/// ------------------------------------------------ ///
	
	
	/// ------------- EVALUATION PARAMETERS ------------ ///
	myFex->storageReliabilityThresholdS2=1250;
	myFex->storageReliabilityThresholdS1=6000;
	
	myFex->numObjectiveOptions = 14;
	myFex->objectiveSpecifier = calloc(myFex->numObjectiveOptions,sizeof(int));
	// myFex->objectiveSpecifierIndex = calloc(myFex->numObjectiveOptions,sizeof(int));
	changeFexSimType(myFex,simulationType);
	myFex->linkCostS2D2 = 0;
	myFex->linkCostS1D2 = 18;
	myFex->linkCostR1S1 = 50;
	myFex->linkCostS1R1 = 0;
	myFex->linkCostGroundBorehole = 40;
	
	/// ------------------------------------------------ ///
	
	/// --------------- POLICY PARAMETERS -------------- ///
	
	//Describe what to normalize inputs between
	myFex->inflowNormMaxMinS2[0] = 200;
	myFex->inflowNormMaxMinS2[1] = 0;
	myFex->inflowNormMaxMinS1[0] = 500;
	myFex->inflowNormMaxMinS1[1] = 0;
	myFex->inflowNormMaxMinR1[0] = 4000;
	myFex->inflowNormMaxMinR1[1] = 0;
	myFex->inflowNormMaxMinR2[0] = 4000;
	myFex->inflowNormMaxMinR2[1] = 0;
	
	myFex->minS2D2 = 12;
	myFex->minS1D2 = 15;
	
	myFex->petNormMaxMinS1[0] = 4.5;
	myFex->petNormMaxMinS1[1] = 0;

	myFex->demandNormMaxMinD2[0] = 73;
	myFex->demandNormMaxMinD2[1] = 36;
	
	//Specify what policy inputs are...
	myFex->numInputOptions = 16;
	//...for cooperative net
	myFex->cooperativeNet = calloc(1,sizeof(struct netOptions));
	myFex->cooperativeNet->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->cooperativeNet->policyInputSpecifier[0] = 1; //S2 storage
	myFex->cooperativeNet->policyInputSpecifier[1] = 1; //S2 inflow
	myFex->cooperativeNet->policyInputSpecifier[2] = 1; //s1 storage
	myFex->cooperativeNet->policyInputSpecifier[3] = 1; //s1 inflow
	myFex->cooperativeNet->policyInputSpecifier[4] = 1; //d2 demand
	myFex->cooperativeNet->policyInputSpecifier[5] = 0; //pet
	myFex->cooperativeNet->policyInputSpecifier[6] = 0; //CM link cap
	myFex->cooperativeNet->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->cooperativeNet->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//
	myFex->cooperativeNet->policyInputSpecifier[9] = 0; //EW link cap
	myFex->cooperativeNet->policyInputSpecifier[10] = 1; //R1 flow
	myFex->cooperativeNet->policyInputSpecifier[11] = 0; //WE link cap
	myFex->cooperativeNet->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->cooperativeNet->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->cooperativeNet->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->cooperativeNet->policyInputSpecifier[15] = 1; //R2 flow
	myFex->cooperativeNet->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	int l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->cooperativeNet->policyInputSpecifier[i] == 1)
		{
			myFex->cooperativeNet->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}

	//...for c2 net
	myFex->c2Net = calloc(1,sizeof(struct netOptions));
	myFex->c2Net->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->c2Net->policyInputSpecifier[0] = 1; //s2 storage
	myFex->c2Net->policyInputSpecifier[1] = 1; //s2 inflow
	myFex->c2Net->policyInputSpecifier[2] = 0; //s1 storage
	myFex->c2Net->policyInputSpecifier[3] = 0; //s1 inflow
	myFex->c2Net->policyInputSpecifier[4] = 1; //d2 demand
	myFex->c2Net->policyInputSpecifier[5] = 0; //pet
	myFex->c2Net->policyInputSpecifier[6] = 0; //CM link cap
	myFex->c2Net->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->c2Net->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//
	myFex->c2Net->policyInputSpecifier[9] = 0; //EW link cap
	myFex->c2Net->policyInputSpecifier[10] = 0; //R1 flow
	myFex->c2Net->policyInputSpecifier[11] = 0; //WE link cap
	myFex->c2Net->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->c2Net->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->c2Net->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->c2Net->policyInputSpecifier[15] = 0; //R2 flow
	myFex->c2Net->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->c2Net->policyInputSpecifier[i] == 1)
		{
			myFex->c2Net->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}

	//...for c1 net
	myFex->c1Net = calloc(1,sizeof(struct netOptions));
	myFex->c1Net->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->c1Net->policyInputSpecifier[0] = 0; //s2 storage
	myFex->c1Net->policyInputSpecifier[1] = 0; //s2 inflow
	myFex->c1Net->policyInputSpecifier[2] = 1; //s1 storage
	myFex->c1Net->policyInputSpecifier[3] = 1; //s1 inflow
	myFex->c1Net->policyInputSpecifier[4] = 0; //d2 demand
	myFex->c1Net->policyInputSpecifier[5] = 0; //pet
	myFex->c1Net->policyInputSpecifier[6] = 0; //CM link cap
	myFex->c1Net->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->c1Net->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//	myFex->c1Net->policyInputSpecifier[9] = 0; //EW link cap
	myFex->c1Net->policyInputSpecifier[10] = 1; //R1 flow
	myFex->c1Net->policyInputSpecifier[11] = 0; //WE link cap
	myFex->c1Net->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->c1Net->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->c1Net->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->c1Net->policyInputSpecifier[15] = 0; //R2 flow
	myFex->c1Net->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->c1Net->policyInputSpecifier[i] == 1)
		{
			myFex->c1Net->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}
	

	//Handy numbers for...
	//...cooperative net
	myFex->cooperativeNet->numNetOutputs = 4;
	myFex->cooperativeNet->numNetInputs = sumInt(myFex->cooperativeNet->policyInputSpecifier,myFex->numInputOptions);
	myFex->cooperativeNet->numLayers = 2 + 1;
	myFex->cooperativeNet->architecture = malloc(myFex->cooperativeNet->numLayers * sizeof(int));
	myFex->cooperativeNet->architecture[0] = myFex->cooperativeNet->numNetInputs;
	myFex->cooperativeNet->architecture[1] = 4;
	// myFex->cooperativeNet->architecture[2] = 3;
	myFex->cooperativeNet->architecture[2] = myFex->cooperativeNet->numNetOutputs;
	myFex->cooperativeNet->numWeights = 0;
	for (i = 0; i < myFex->cooperativeNet->numLayers - 1; i++)
	{
		myFex->cooperativeNet->numWeights = myFex->cooperativeNet->numWeights + myFex->cooperativeNet->architecture[i]*myFex->cooperativeNet->architecture[i+1] + myFex->cooperativeNet->architecture[i+1];
	}
	myFex->cooperativeNet->betaInd = 0;
	for (j = 0; j < (myFex->cooperativeNet->numLayers - 2); j++)
	{
		myFex->cooperativeNet->betaInd = myFex->cooperativeNet->betaInd + myFex->cooperativeNet->architecture[j]*myFex->cooperativeNet->architecture[j+1];
	}
	myFex->cooperativeNet->betaIndMax = myFex->cooperativeNet->betaInd;
	for (j = 1; j < (myFex->cooperativeNet->numLayers - 1);j++)
	{
		myFex->cooperativeNet->betaIndMax = myFex->cooperativeNet->betaIndMax + myFex->cooperativeNet->architecture[j];
	}
	
	/// ------------------------------------------------ ///
	
	/// ---------------- ALLOCATE MEMORY --------------- /// -- would be good to also have a reallocate function (i.e. frees memory, then allocates to a new fexT
	myFex->warmup = 10000;
	/*myFex->fexFlowsS2 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsS1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsR1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsR2 = calloc(myFex->fexT,sizeof(double));
	myFex->s1FishReleases = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS2D2TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS1D2TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapR1S1TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS1R1TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapBoreholeTimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->fexPets= calloc(myFex->fexT,sizeof(double));
	// myFex->fexDemandsc1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexDemandsD2 = calloc(myFex->fexT,sizeof(double));
	myFex->fexDoys = calloc(myFex->fexT,sizeof(int));*/
	/// ------------------------------------------------ ///
}

void matXvec(double *mat, double *vec, int M, int N, double *product)
{
	//Performs matrix multipilcation of a [MxN] matrix (flattened such that mat[i][j] is equivalent to mat[j*M + i]) with a [N] vector (vec) and stores result in pointer product (a [M] vector)
	int i,j;
	for (i=0;i<M;i++)
	{
		product[i] = 0;
		for (j=0;j<N;j++)
		{
			product[i] = product[i] + mat[j*M + i]*vec[j];
		}
	}
}

int poisRnd(int lambda)
{
	int n = 0;
	double limit = exp(-lambda);
	double x = (double)rand()/(double)RAND_MAX;
	while (x > limit)
	{
		n++;
		x *= (double)rand()/(double)RAND_MAX;
	}
	return n;
}


double normRnd (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}
FILE * readParm(FILE *f, struct parma *myParma, int *offset)
{
	int i;
	char line[256];
	fgets(line,sizeof(line),f);
	myParma->numHarmonics = atoi(line);
	
	
	myParma->harmonics = malloc( myParma->numHarmonics * sizeof(double));
	for (i = 0; i < myParma->numHarmonics; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->harmonics[i] = atof(line);
		
	}
	
	myParma->thetas = malloc( (2*myParma->numHarmonics + 1) * sizeof(double));
	for (i = 0; i < (2*myParma->numHarmonics + 1); i++)
	{
		fgets(line,sizeof(line),f);
		myParma->thetas[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->numAR = atoi(line);
	*offset = fmax(myParma->numAR,*offset);
	
	myParma->AR = malloc( myParma->numAR * sizeof(double));
	for (i = 0; i < myParma->numAR ; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->AR[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->numMA = atoi(line);
	*offset = fmax(myParma->numMA,*offset);
	
	myParma->MA = malloc( myParma->numMA * sizeof(double));
	for (i = 0; i < myParma->numMA ; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->MA[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->cons = atof(line);
	
	fgets(line,sizeof(line),f);
	myParma->varc = atof(line);
	fgets(line,sizeof(line),f);
	if (line[0] != '#')
	{
		printf("warning: formatting problem in parma file\n");
	}
	
	return f;
}
FILE * readTheta(FILE *f, struct parma *myParma)
{
	int i;
	char line[256];
	fgets(line,sizeof(line),f);
	myParma->numHarmonics = atoi(line);
	
	
	myParma->harmonics = malloc( myParma->numHarmonics * sizeof(double));
	for (i = 0; i < myParma->numHarmonics; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->harmonics[i] = atof(line);
		
	}
	
	myParma->thetas = malloc( (2*myParma->numHarmonics + 1) * sizeof(double));
	for (i = 0; i < (2*myParma->numHarmonics + 1); i++)
	{
		fgets(line,sizeof(line),f);
		myParma->thetas[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	if (line[0] != '#')
	{
		printf("warning: formatting problem in theta file\n");
	}
	
	return f;
}

void createEnsemble(struct fexEnsemble *myEnsemble, int framing)
{
	//AN ENSEMBLE (fexEnsemble) CONTAINS EVERYTHING NOT USED IN THE SIMULATION AND POINTERS (fexEnsemble.ensemble) TO THE PROBLEMS
	//A PROBLEM (fexProblem) CONTAINS ONLY WHAT IS NEEDED FOR THE SIMULATION
	
	int i, j;
	myEnsemble->NUM_AUTOCORRELATED_PROCESSES = 6;
	

	char parmaFid[256];
	char choleskyFid[256];
	if (framing == 1 )
	{
		sprintf(parmaFid,"parma00.arma");
		sprintf(choleskyFid,"cholesky00.arma");
	}
	else
	{
		sprintf(parmaFid,"parma20.arma");
		sprintf(choleskyFid,"cholesky20.arma");
	}
	
	//edit: 18/05/18 - Add variable storage threshold
	char s2ThetaFid[256];
	char s1ThetaFid[256];
	
	sprintf(s2ThetaFid,"s2_storage.thetas");
	sprintf(s1ThetaFid,"s1_storage.thetas");
	//endedit

	/// --------------- READ CHOLESKY ------------------ ///
	myEnsemble->cholesky = malloc( myEnsemble->NUM_AUTOCORRELATED_PROCESSES * (sizeof(double*)) );
	FILE *fChol = fopen(choleskyFid,"r");
	for (i=0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		myEnsemble->cholesky[i] = malloc( myEnsemble->NUM_AUTOCORRELATED_PROCESSES * (sizeof(double)) );
		for (j=0; j < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; j++)
		{
			if (!fscanf(fChol, "%lf", &myEnsemble->cholesky[i][j]))
			break;
		}
	}
	fclose(fChol);
	/// ----------- READ AND ALLOCATE PARMA ------------ ///
	
	FILE *fParm = fopen(parmaFid,"r");
	myEnsemble->offset = 0;
	fParm = readParm(fParm,&myEnsemble->S2Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->S1Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->R1Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->R2Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->petParma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->D2Parma,&myEnsemble->offset);
	fclose(fParm);

	//edit: 18/05/18 - Add variable storage threshold
	FILE *fS2 = fopen(s2ThetaFid,"r");
	fS2 = readTheta(fS2,&myEnsemble->s2Theta);
	fclose(fS2);
	
	FILE *fS1 = fopen(s1ThetaFid,"r");
	fS1 = readTheta(fS1,&myEnsemble->s1Theta);
	fclose(fS1);
	//endedit

	/// --------------- BREAK PARAMETERS --------------- ///
	if (framing == 4)
	{
		myEnsemble->breakLambdaProbabilityS2D2[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS2D2[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1D2[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1D2[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityR1S1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityR1S1[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1R1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1R1[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityBorehole[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityBorehole[1] = 1000000000000; //lambda for poisson duration between break
	}
	else
	{
		myEnsemble->breakLambdaProbabilityS2D2[0] = 5; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS2D2[1] = 800; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1D2[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1D2[1] = 300; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityR1S1[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityR1S1[1] = 300; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1R1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1R1[1] = 0; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityBorehole[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityBorehole[1] = 500; //lambda for poisson duration between break
	}
	/// ------------------------------------------------ ///

	
	/// ------------- S1 FISHERIES ------------- ///
	
	//s1FisheriesAnnual will be uniformly released over a random length between min/maxDailyFishRelease during start/endAllowableFishReleaseDoy
	myEnsemble->s1FisheriesAnnual = 900; //Ml/y
	if (framing == 4)
	{
		myEnsemble->startAllowableFishReleaseDoy = 60; //i.e. March 1
		myEnsemble->endAllowableFishReleaseDoy = 365; //i.e. Dec 31
	}
	else
	{
		myEnsemble->startAllowableFishReleaseDoy = 244; //i.e. Sept 1
		myEnsemble->endAllowableFishReleaseDoy = 273; //i.e. Sept 31	
	}
	myEnsemble->maxFishReleaseLength = 14; //Entire annual allowance must be released over a maximum of 14 consecutive days
	myEnsemble->maxDailyFishRelease = 360; //Ml/d - note- that since this isn't divisible by 900, the actual max is 300
	myEnsemble->minFishReleaseLength = (int)ceil(myEnsemble->s1FisheriesAnnual/myEnsemble->maxDailyFishRelease);
	
	/// ------------------------------------------------ ///
	
	/// ---------------- CREATE ENSEMBLE --------------- ///
	myEnsemble->ensemble = calloc(myEnsemble->ensembleSize,sizeof(struct fexProblem*));
	srand(myEnsemble->baseSeed);
	for (i=0;i<myEnsemble->ensembleSize;i++)
	{
		myEnsemble->ensemble[i] = calloc(1,sizeof(struct fexProblem));
		createFex(myEnsemble->ensemble[i], myEnsemble, myEnsemble->ensembleT, myEnsemble->baseSeed,myEnsemble->simulationType);
		int seed = rand()+i;
		myEnsemble->ensemble[i]->parent = myEnsemble;
		myEnsemble->ensemble[i]->seed = seed;
	}
	/// ------------------------------------------------ ///

	
	/// ---------------- CHANGE FRAMING --------------- ///
	myEnsemble->framing = framing;

}

void destroyFex(struct fexProblem *myFex)
{
	
	//free(myFex->fexFlowsS2);
	//free(myFex->fexFlowsS1);
	//free(myFex->fexFlowsR1);
	//free(myFex->fexFlowsR2);
	//
	//free(myFex->s1FishReleases);
	//
	//free(myFex->fexDoys);
	//
	//free(myFex->linkCapS2D2TimeVarying);
	//free(myFex->linkCapS1D2TimeVarying);
	//free(myFex->linkCapR1S1TimeVarying);
	//free(myFex->linkCapS1R1TimeVarying);
	//free(myFex->linkCapBoreholeTimeVarying);
	//
	//free(myFex->fexPets);	
	//free(myFex->fexDemandsD2);
	// free(myFex->fexDemandsc1);
	
	free(myFex->objectiveSpecifier);
	
	free(myFex->cooperativeNet->policyInputSpecifier);
	free(myFex->cooperativeNet->policyInputSpecifierIndex);
	free(myFex->cooperativeNet->architecture);
	free(myFex->cooperativeNet);	
	free(myFex->c2Net->policyInputSpecifier);
	free(myFex->c2Net->policyInputSpecifierIndex);
	free(myFex->c2Net->architecture);
	free(myFex->c2Net);
	free(myFex->c1Net->policyInputSpecifier);
	free(myFex->c1Net->policyInputSpecifierIndex);
	free(myFex->c1Net->architecture);
	free(myFex->c1Net);
}

void destroyParma(struct parma *myParma)
{
	free(myParma->harmonics);
	free(myParma->thetas);
	free(myParma->AR);
	free(myParma->MA);
}
void destroyEnsemble(struct fexEnsemble *myEnsemble)
{
	//Free fexEnsemble (and fex problems within fexEnsemble.ensemble)
	int i,j;
	for (i = 0;i<myEnsemble->ensembleSize;i++)
	{
		destroyFex(myEnsemble->ensemble[i]);
		for (j=0;j<myEnsemble->ensembleSize;j++)
		{
			if (myEnsemble->ensemble[i]->parent == myEnsemble->ensemble[j]->parent && i != j)
			{
				myEnsemble->ensemble[i]->parent = NULL;
				break;
			}
		}
		
		if (myEnsemble != myEnsemble->ensemble[i]->parent && myEnsemble->ensemble[i]->parent != NULL)
		{
			
			destroyEnsemble(myEnsemble->ensemble[i]->parent);
			free(myEnsemble->ensemble[i]->parent);
		}
	}
	for (i=0 ; i< myEnsemble->ensembleSize;i++)
	{
		free(myEnsemble->ensemble[i]);
	}
	
	free(myEnsemble->ensemble);
	
	destroyParma(&myEnsemble->S2Parma);
	destroyParma(&myEnsemble->S1Parma);
	destroyParma(&myEnsemble->R2Parma);
	destroyParma(&myEnsemble->R1Parma);
	destroyParma(&myEnsemble->D2Parma);
	destroyParma(&myEnsemble->petParma);
	
	for (i = 0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		free(myEnsemble->cholesky[i]);
	}
	free(myEnsemble->cholesky);
}
void evaluateNet(double *outputs, double *inputs, int *architecture, int numLayers, double *weights)
{
	//evaluate an RBF net with described architecture and inputs, storing result in outputs
	//TOTAL NUMBER OF WEIGHTS = (numNetInputs+1)*numNetHidden + numNetOutputs*(numNetHidden+1)
	//Implements rbfn according to http://mccormickml.com/2013/08/15/radial-basis-function-network-rbfn-tutorial/

	int i,j,k,m;
	double *phis_old_ptr;
	for (k = 1; k < numLayers - 1; k++)
	{
		// double *phis_old;
		if (k==1)
		{
			phis_old_ptr = inputs;
		}
		// phis_old = malloc(architecture[k-1] * sizeof(double));
		// for (i = 0; i < architecture[k-1] ; i++)
		// {
			// phis_old[i] = phis_old_ptr[i];
		// }
		
		double *phis = malloc(architecture[k] * sizeof(double));
		for (i = 0; i < architecture[k]; i++)
		{
			double sqrD = 0;
			for (j = 0; j < architecture[k-1]; j++)
			{
				int ind = j + i*architecture[k-1];
				for (m = 0; m < (k-1); m++)
				{
					ind = ind + architecture[m]*architecture[m+1];
				}
				// printf("ind: %d, weight: %f\n",ind,weights[ind]);
				sqrD = sqrD + weights[ind]*phis_old_ptr[j];
			}
			int ind = i;
			for (j = 1; j < k; j++)
			{
				ind = ind + architecture[j];
			}
			for (j = 0; j < (numLayers - 2); j++)
			{
				ind = ind + architecture[j]*architecture[j+1];
			}
			// printf("ind: %d, weight: %f\n",ind,weights[ind]);
			phis[i] = weights[ind] + 2/(1+exp(-2*sqrD)) - 1;;
		}
		if (k != 1)
		{
			free(phis_old_ptr);
		}
		phis_old_ptr = phis;
	}
	for (i=0; i < architecture[numLayers - 1]; i++)
	{
		int ind = i;
		for (j = 1 ; j < numLayers - 1; j++)
		{
			ind = ind + architecture[j] + architecture[j-1]*architecture[j];
		}
		double Y = weights[ind];
		// printf("ind: %d, weight: %f\n",ind,weights[ind]);
		for (j=0; j < architecture[numLayers-2]; j++)
		{
			Y = Y + phis_old_ptr[j]*weights[ind + architecture[numLayers-1] + j + i*(architecture[numLayers-2]-1)];
			// printf("ind: %d, weight: %f\n",ind + architecture[numLayers-1] + j + i*(architecture[numLayers-2]-1),weights[ind + architecture[numLayers-1] + j + i*(architecture[numLayers-2]-1)]);
		}
	
		if (Y < -1)
		{
			outputs[i] = 0;
		}
		else if (Y < 1)
		{
			outputs[i] = (Y + 1)/2;
		}
		else
		{
			outputs[i] = 1;
		}
	}
	free(phis_old_ptr);
	
	fflush(stdout);
}
double sumDouble(double *vec, int T)
{
	//Sum of vector (vec) of doubles (length = T)
	int i;
	double sum = 0;
	for (i=0;i<T;i++)
	{
		sum = sum + vec[i];
	}
	return sum;
}

//edit: loram - 4/Apr/18
double generateArmaOneStep(struct parma *myParma, double *laggedResidual, double *laggedInnovation, double innovation, int T, int offset)
{
	int j;
	double ar = 0;
	double sigma = sqrt(myParma->varc);
	for (j = 0; j < myParma->numAR; j++)
	{
		ar = ar + laggedResidual[offset - myParma->numAR + j]*myParma->AR[j]; //assuming here that AR[0] is the parameter for the residual myParma->numAR time-steps ago
	}
	double ma = 0;
	for (j = 0; j < myParma->numMA; j++)
	{
		ma = ma + laggedInnovation[offset - myParma->numMA+ j]*myParma->MA[j];

	}
	double residual = myParma->cons + ar + ma + innovation*sigma;
	return residual;
}
void randomCorrelatedNormal(int N, double **cholesky, double *randomCorrNormal)
{
	int i,k;
	double *randomNormal = malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		randomNormal[i] = normRnd(0,1);
	}
	for (i = 0; i < N; i++)
	{
		randomCorrNormal[i] = 0;
		for (k = 0; k < N; k++)
		{
			randomCorrNormal[i] = randomCorrNormal[i] + randomNormal[k]*cholesky[k][i];
		}
	}
	free(randomNormal);
}


void simFex(double* weights, double* objs, double* consts,struct fexProblem *myFex)
{
	// simulate a fexProblem for weights, storing objectives in objs
	// consts is currently unused
	
	/// ---------------- ALLOCATE MEMORY --------------- ///
	
	int i, j, t, l;
	double storageS2 = 0; //S2 stor
	double spillS2 = 0; //S2 spill
	double storageS1 = 0; //s1 stor
	double spillS1 = 0; //s1 spill
	double supplyToD2 = 0; //How much to supply to D2
	double supplyToC1 = 0; //How much to supply to D1
	double volumeS2D2 = 0; //s2 release
	double volumeS1D2 = 0; //s1 release
	double volumeR1S1 = 0; //river R1 pumped abstraction to s1
	double volumeS1R1 = 0; //s1all release to R1 for DE2
	double volumeBorhole = 0; //c1 groundwater supplement
	
	double fexFlowS1;
	double fexFlowS2;
	double fexFlowR1;
	double fexFlowR2;
	double fexPet;
	double fexDemandD1;
	double fexDemandD2;
	double s1FishReleases;

	double elevation; double area; double evap;
	double S2D2AnnualTotal = 0;
	double S1D2AnnualTotal = 0;
	double R1S1AnnualTotal = 0;
	double S2D2MonthlyTotal = 0;
	double S1R1AnnualTotal = 0;
	/*
	double *pumpCostTotal = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *pumpCostC2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *pumpCostC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVc2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVtotal = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitSqrD2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitLinD2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitSqrC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitLinC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *gsSqrC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *gsLinC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *storageReliabilityS2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *storageReliabilityS1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	*/
	double gsLinC1 = 0;
	double gsSqrC1 = 0;
	double deficitSqrC1 = 0;
	double deficitLinC1 = 0;
	double deficitSqrD2 = 0;
	double deficitLinD2 = 0;
	double storageReliabilityS2 = 0;
	double storageReliabilityS1 = 0;
	double NPVc2 = 0;
	double NPVC1 = 0;
	double NPVtotal = 0;
	double pumpCostC2 = 0;
	double pumpCostC1 = 0;
	double pumpCostTotal = 0;

	const double DOUBLE_PRECISION = pow(2,-52);
	const int DAYS_IN_YEAR = 365;
	const int DAYS_IN_MONTH = 30;
	const double MM_M2_TO_ML = pow(10,-6);
	// NOTE: outputting objectives over time through use of 'timeObjectives' and 'includetimeObjectives' is currently untested and commented out
	// timeObjectives = calloc(myFex->warmup + myFex->numObjectives,sizeof(double)); // NOTE: if includetimeObjectives == 1, then this data and the associated objectives will not be freed (not even in destroyFex)
	
	/// ------------------------------------------------ ///
	
	
	/// ---------------- SETUP PRINTING ---------------- ///
	FILE *fid;
	if (myFex->printDynamics == 1)
	{
		fid = fopen(myFex->dynamicFid,"w");
		fprintf(fid,"//Sc Sw Ic Iw Iex Ith Dm Ds ucm uwm uew uwe gs wc ww doy evapC evapW s2Env s1Env");
//		for (i = 0; i < myFex->numObjectives; i++)
//		{
//			fprintf(fid, " obj%d", i);
//		}
		fprintf(fid,"\n");
	}
	/// ------------------------------------------------ ///
	
	
	/// ---------------- SETUP NETWORKS ---------------- ///
	
	double *c2Weights;
	double *C1Weights;
	int numNetInputs;
	int numNetOutputs;
	int *architecture;
	int numLayers;
	int numWeights;
	int *policyInputSpecifier;
	int *policyInputSpecifierIndex;
	switch (myFex->simulationType)
	{
		case 0:
			numNetInputs = myFex->cooperativeNet->numNetInputs;
			numNetOutputs = myFex->cooperativeNet->numNetOutputs;
			numLayers = myFex->cooperativeNet->numLayers;
			architecture = myFex->cooperativeNet->architecture;
			numWeights = myFex->cooperativeNet->numWeights;
			policyInputSpecifier = myFex->cooperativeNet->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->cooperativeNet->policyInputSpecifierIndex;
			break;
		case 1:
			numNetInputs = myFex->c2Net->numNetInputs;
			numNetOutputs = myFex->c2Net->numNetOutputs;
			numWeights = myFex->c2Net->numWeights;
			policyInputSpecifier = myFex->c2Net->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->c2Net->policyInputSpecifierIndex;
			break;
		case 2:
			numNetInputs = myFex->c1Net->numNetInputs;
			numNetOutputs = myFex->c1Net->numNetOutputs;
			numWeights = myFex->c1Net->numWeights;
			policyInputSpecifier = myFex->c1Net->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->c1Net->policyInputSpecifierIndex;
			break;
		case 3:
			c2Weights = calloc(myFex->c2Net->numWeights,sizeof(double));
			for (i = 0; i < myFex->c2Net->numWeights;i++)
			{
				c2Weights[i] = weights[i];
			}
			C1Weights = calloc(myFex->c1Net->numWeights,sizeof(double));
			for (j = 0; j< myFex->c1Net->numWeights;j++)
			{
				C1Weights[j] = weights[i+j];
			}	
			numNetOutputs = myFex->cooperativeNet->numNetOutputs;
			break;
	}		
		//edit: 18/05/18 - Add variable storage threshold
	double reliabilityMod = 1;
	if (myFex->parent->framing == 3)
	{
		reliabilityMod = 0.5;
	}
	//endedit
	/// ------------------------------------------------ ///
	
	//edit: LoRam 4/Apr/18
	//Initialize simulation
	// collect parmas in correct order
	struct parma *allParma = calloc(myFex->parent->NUM_AUTOCORRELATED_PROCESSES, sizeof(struct parma));
	allParma[0] = myFex->parent->S2Parma;
	allParma[1] = myFex->parent->S1Parma;
	allParma[2] = myFex->parent->R1Parma;
	allParma[3] = myFex->parent->R2Parma;
	allParma[4] = myFex->parent->petParma;
	allParma[5] = myFex->parent->D2Parma;
	
	//warmup
	srand(myFex->seed);
	double **innovationHolders = calloc(myFex->parent->NUM_AUTOCORRELATED_PROCESSES, sizeof(double*));
	double **residualHolders = calloc(myFex->parent->NUM_AUTOCORRELATED_PROCESSES, sizeof(double*));
	for (i = 0; i < myFex->parent->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		innovationHolders[i] = calloc(myFex->parent->offset,sizeof(double));
		residualHolders[i] = calloc(myFex->parent->offset,sizeof(double));
		for (j = 0; j < myFex->parent->offset; j++)
		{
			residualHolders[i][j] = 0;
			innovationHolders[i][j] = 0;
		}
	}
	
	for (t = 0; t < myFex->warmup; t++)
	{
		double *randomCorrNormal = calloc(myFex->parent->NUM_AUTOCORRELATED_PROCESSES,sizeof(double));
		randomCorrelatedNormal(myFex->parent->NUM_AUTOCORRELATED_PROCESSES, myFex->parent->cholesky, randomCorrNormal);
		for (i = 0; i < myFex->parent->NUM_AUTOCORRELATED_PROCESSES; i++)
		{
			double residual = generateArmaOneStep(&allParma[i], residualHolders[i], innovationHolders[i], randomCorrNormal[i], 1, myFex->parent->offset);
			for (j = 1; j < myFex->parent->offset; j++)
			{
				residualHolders[i][j-1] = residualHolders[i][j];
				innovationHolders[i][j-1] = innovationHolders[i][j];
			}
			residualHolders[i][j-1] = residual;
			innovationHolders[i][j-1] = randomCorrNormal[i];
		}
		free(randomCorrNormal);
	}
	
	storageS2 = ((double)rand()/(double)RAND_MAX)*(myFex->storCapS2 - myFex->deadStorS2) + myFex->deadStorS2;
	storageS1 = ((double)rand()/(double)RAND_MAX)*(myFex->storCapS1 - myFex->deadStorS1) + myFex->deadStorS1;
	
	int startOfFishReleases = myFex->fexT;
	int lengthOfFishReleases = 0;
	double averageFishReleases = 0;
	
	int NUM_POSSIBLE_PIPE_FAILURES = 5;
	int *pipeCondition = calloc(NUM_POSSIBLE_PIPE_FAILURES,sizeof(int));
	int *timeUntilBreak = calloc(NUM_POSSIBLE_PIPE_FAILURES,sizeof(int));
	double *breakLambdas = calloc(NUM_POSSIBLE_PIPE_FAILURES * 2,sizeof(double));
	for (i = 0; i < NUM_POSSIBLE_PIPE_FAILURES; i++)
	{
		pipeCondition[i] = 0;
		switch (i)
		{
			case 0:
				breakLambdas[i*2] = myFex->parent->breakLambdaProbabilityS2D2[0];
				breakLambdas[i*2 + 1] = myFex->parent->breakLambdaProbabilityS2D2[1];
			case 1:
				breakLambdas[i*2] = myFex->parent->breakLambdaProbabilityS1D2[0];
				breakLambdas[i*2 + 1] = myFex->parent->breakLambdaProbabilityS1D2[1];
			case 2:
				breakLambdas[i*2] = myFex->parent->breakLambdaProbabilityR1S1[0];
				breakLambdas[i*2 + 1] = myFex->parent->breakLambdaProbabilityR1S1[1];
			case 3:
				breakLambdas[i*2] = myFex->parent->breakLambdaProbabilityS1R1[0];
				breakLambdas[i*2 + 1] = myFex->parent->breakLambdaProbabilityS1R1[1];
			case 4:
				breakLambdas[i*2] = myFex->parent->breakLambdaProbabilityBorehole[0];
				breakLambdas[i*2 + 1] = myFex->parent->breakLambdaProbabilityBorehole[1];
		}
		timeUntilBreak[i] = poisRnd(breakLambdas[i*2+1]);
	}
	
	int doy = (rand()+1)%365;
	if (doy == 0)
	{
		doy = 365;
	}
	//endedit

	/// ------------------------------------------------ ///
	
	//Begin simulation
	for (t=0;t<(myFex->warmup+ myFex->fexT);t++)
	{
		//edit: LoRam 4/Apr/18
		doy++;
		if (doy == 366)
		{
			doy = 1;
		}
		//Generate Forcing
		double *randomCorrNormal = calloc(myFex->parent->NUM_AUTOCORRELATED_PROCESSES,sizeof(double));
		randomCorrelatedNormal(myFex->parent->NUM_AUTOCORRELATED_PROCESSES, myFex->parent->cholesky, randomCorrNormal);
		for (i = 0; i < myFex->parent->NUM_AUTOCORRELATED_PROCESSES; i++)
		{
			double residual = generateArmaOneStep(&allParma[i], residualHolders[i], innovationHolders[i], randomCorrNormal[i], 1, myFex->parent->offset);
			for (j = 1; j < myFex->parent->offset; j++)
			{
				residualHolders[i][j-1] = residualHolders[i][j];
				innovationHolders[i][j-1] = innovationHolders[i][j];
			}
			residualHolders[i][j-1] = residual;
			innovationHolders[i][j-1] = randomCorrNormal[i];
			
			double periodic = allParma[i].thetas[0];
			for (j = 0; j < allParma[i].numHarmonics; j++)
			{
				periodic = periodic + allParma[i].thetas[(j+1)*2-1]*sin(allParma[i].harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
				periodic = periodic + allParma[i].thetas[(j+1)*2]*cos(allParma[i].harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
			}
			switch (i)
			{
				case 0:
					fexFlowS2 = fmax(periodic*exp(residual),0);
				case 1:
					fexFlowS1 = fmax(periodic*exp(residual),0);
				case 2:
					fexFlowR1 = fmax(periodic*exp(residual),0);
				case 3:
					fexFlowR2 = fmax(periodic*exp(residual),0);
				case 4:
					fexPet = fmax(periodic + residual,0);
				case 5:
					fexDemandD2 = fmax(periodic + residual,0);
			}
				
		}
		free(randomCorrNormal);
		//edit: 18/05/18 - Add variable storage threshold
		double storageThresholdS2 = myFex->parent->s2Theta.thetas[0];
		for (j = 0; j < myFex->parent->s2Theta.numHarmonics;j++)
		{
			storageThresholdS2 = storageThresholdS2 + myFex->parent->s2Theta.thetas[(j+1)*2-1]*sin(myFex->parent->s2Theta.harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
			storageThresholdS2 = storageThresholdS2 + myFex->parent->s2Theta.thetas[(j+1)*2]*cos(myFex->parent->s2Theta.harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
		}
		double storageThresholdS1 = myFex->parent->s1Theta.thetas[0];
		for (j = 0; j < myFex->parent->s1Theta.numHarmonics;j++)
		{
			storageThresholdS1 = storageThresholdS1 + myFex->parent->s1Theta.thetas[(j+1)*2-1]*sin(myFex->parent->s1Theta.harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
			storageThresholdS1 = storageThresholdS1 + myFex->parent->s1Theta.thetas[(j+1)*2]*cos(myFex->parent->s1Theta.harmonics[j]*M_PI*doy/DAYS_IN_YEAR);
		}
		//endedit
		
		if (doy == myFex->parent->startAllowableFishReleaseDoy)
		{
			lengthOfFishReleases = rand()%(myFex->parent->maxFishReleaseLength - myFex->parent->minFishReleaseLength + 1);
			lengthOfFishReleases = lengthOfFishReleases + myFex->parent->minFishReleaseLength;
			averageFishReleases = myFex->parent->s1FisheriesAnnual/lengthOfFishReleases;
			startOfFishReleases = rand()%(myFex->parent->endAllowableFishReleaseDoy - myFex->parent->startAllowableFishReleaseDoy - lengthOfFishReleases + 2);
		}
		if (startOfFishReleases < 1)
		{
			if (lengthOfFishReleases == 0)
			{
				startOfFishReleases = myFex->fexT;
				s1FishReleases = 0;
			}
			else
			{
				s1FishReleases = averageFishReleases;
				lengthOfFishReleases--;
			}
		}
		else
		{
			s1FishReleases = 0;
		}
		startOfFishReleases--;
		
		// breaks
		for (i = 0; i < NUM_POSSIBLE_PIPE_FAILURES; i++)
		{
			if (pipeCondition[i] > 0)
			{
				pipeCondition[i]--;
				if (pipeCondition[i] == 0)
				{
					timeUntilBreak[i] = poisRnd(breakLambdas[i*2+1]);
				}
			}
			if (timeUntilBreak[i] < 1)
			{
				pipeCondition[i] = poisRnd(breakLambdas[i*2]);
			}
			timeUntilBreak[i]--; 
		}
		//
		//endedit		

		/// ---------- APPLY PHYSICAL FORCING ---------- ///
		double s2Env;
		double s1Env;

		storageS2 = storageS2 + fexFlowS2 - myFex->s2CompensationFixed;
		storageS1 = storageS1 + fexFlowS1 - myFex->s1CompensationFixed - s1FishReleases;		
		
				if (storageS2 < 0)
		{
			s2Env = fmax(myFex->s2CompensationFixed + storageS2,0);
			storageS2 = 0;
		}
		else
		{
			s2Env = myFex->s2CompensationFixed;
		}
		if (storageS1 < 0)
		{
			s1Env = fmax(myFex->s1CompensationFixed + s1FishReleases + storageS1,0);
			storageS1 = 0;
		}
		else
		{
			s1Env = myFex->s1CompensationFixed + s1FishReleases;
		}


		elevation = myFex->storElevationS2[0]*(pow(storageS2,myFex->storElevationS2[1]));
		area = myFex->elevationAreaS2[0]*(pow(elevation,myFex->elevationAreaS2[1]));
		
		evap = fexPet;

		double evapS2 = (evap*area)*MM_M2_TO_ML;
		storageS2 = fmax((storageS2 - evapS2),0);
		
		elevation = myFex->storElevationS1[0]*(pow(storageS1,myFex->storElevationS1[1]));
		area = myFex->elevationAreaS1[0]*(pow(elevation,myFex->elevationAreaS1[1]));
		evap = fexPet;
		double evapS1 = (evap*area)*MM_M2_TO_ML;
		storageS1 = fmax((storageS1 - evapS1),0);
	
		/// -------------------------------------------- ///
		
		
		/// ------- CREATE AND EVALUATE POLICIES ------- ///
		
		double *targetReleases = calloc(numNetOutputs,sizeof(double));
		if (myFex->simulationType != 3)
		{	
			double *netInputs = calloc(numNetInputs,sizeof(double));
			if (policyInputSpecifier[0] == 1)
			{			
				//Normalize storage
				netInputs[policyInputSpecifierIndex[0]] = (storageS2 - myFex->deadStorS2)/(myFex->storCapS2 - myFex->deadStorS2 + myFex->inflowNormMaxMinS2[0]);
			}
			if (policyInputSpecifier[1] == 1)
			{
				//Normalize s2 inflows
				netInputs[policyInputSpecifierIndex[1]] = (fexFlowS2 - myFex->inflowNormMaxMinS2[1])/(myFex->inflowNormMaxMinS2[0] - myFex->inflowNormMaxMinS2[1]);
			}
			
			if (policyInputSpecifier[2] == 1)
			{			
				//Normalize s1 storage
				netInputs[policyInputSpecifierIndex[2]] = (storageS1 - myFex->deadStorS1)/(myFex->storCapS1 - myFex->deadStorS1 + myFex->inflowNormMaxMinS1[0]);
			}
			if (policyInputSpecifier[3] == 1)
			{
				//Normalize s1 inflows
				netInputs[policyInputSpecifierIndex[3]] = (fexFlowS1 - myFex->inflowNormMaxMinS1[1])/(myFex->inflowNormMaxMinS1[0] - myFex->inflowNormMaxMinS1[1]);
			}
			
			if (policyInputSpecifier[4] == 1)
			{
				//Normalize d2down demand
				netInputs[policyInputSpecifierIndex[4]] = (fexDemandD2 - myFex->demandNormMaxMinD2[1])/(myFex->demandNormMaxMinD2[0] - myFex->demandNormMaxMinD2[1]);
			}
			if (policyInputSpecifier[5] == 1)
			{
				//Normalize pet
				netInputs[policyInputSpecifierIndex[5]] = (fexPet - myFex->petNormMaxMinS1[1])/(myFex->petNormMaxMinS1[0] - myFex->petNormMaxMinS1[1]);
			}
			if (policyInputSpecifier[6] == 1)
			{
				if (pipeCondition[0] > 0)
				{
					netInputs[policyInputSpecifierIndex[6]] = 0;
				}
				else
				{
					netInputs[policyInputSpecifierIndex[6]] = 1;
				}
			}
			if (policyInputSpecifier[7] == 1)
			{
				if (pipeCondition[1] > 0)
				{
					netInputs[policyInputSpecifierIndex[7]] = 0;
				}
				else
				{
					netInputs[policyInputSpecifierIndex[7]] = 1;
				}
			}			
			if (policyInputSpecifier[8] == 1)
			{
				//Normalize C1 demand
			}
			if (policyInputSpecifier[9] == 1)
			{
				//Normalize r1 release capacity
				if (pipeCondition[2] > 0)
				{
					netInputs[policyInputSpecifierIndex[9]] = 0;
				}
				else
				{
					netInputs[policyInputSpecifierIndex[9]] = 1;
				}
			}
			if (policyInputSpecifier[10] == 1)
			{
				//Normalize R1 inflows
				netInputs[policyInputSpecifierIndex[10]] = (fexFlowR1 - myFex->inflowNormMaxMinR1[1])/(myFex->inflowNormMaxMinR1[0] - myFex->inflowNormMaxMinR1[1]);
			}
			if (policyInputSpecifier[11] == 1)
			{
				//Normalize s1 release to R1 capacity
				if (pipeCondition[3] > 0)
				{
					netInputs[policyInputSpecifierIndex[11]] = 0;
				}
				else
				{
					netInputs[policyInputSpecifierIndex[11]] = 1;
				}
			}
			if (policyInputSpecifier[12] == 1)
			{
				if (pipeCondition[4] > 0)
				{
					netInputs[policyInputSpecifierIndex[12]] = 0;
				}
				else
				{
					netInputs[policyInputSpecifierIndex[12]] = 1;
				}
			}
			if (policyInputSpecifier[13] == 1)
			{
				//time component sin
				netInputs[policyInputSpecifierIndex[13]] = sin(2*M_PI*doy/365);
			}
			if (policyInputSpecifier[14] == 1)
			{
				//time component cos
				netInputs[policyInputSpecifierIndex[14]] = cos(2*M_PI*doy/365);
			}
			if (policyInputSpecifier[15] == 1)
			{
				//normalize R2 flow
				netInputs[policyInputSpecifierIndex[15]] = (fexFlowR2 - myFex->inflowNormMaxMinR2[1])/(myFex->inflowNormMaxMinR2[0] - myFex->inflowNormMaxMinR2[1]);
			}
			for (i = 0; i < numNetInputs;i++)
			{
				netInputs[i] = fmax(0, fmin(netInputs[i],1)); //In case variables have exceeded their normalization limits
			}
			evaluateNet(targetReleases,netInputs,architecture,numLayers,weights);
			
			free(netInputs);
		}
		double targS2D2, targS1D2, targR1S1, targS1R1;
		if (myFex->simulationType != 2)
		{ 
			targS2D2 = targetReleases[0];
			targS1D2 = targetReleases[1];
			if (myFex->simulationType != 1)
			{
				targR1S1 = targetReleases[2];
				targS1R1 = targetReleases[3];
			}
			else 
			{
				targR1S1 = 0;
				targS1R1 = 0;
			}
		}
		else
		{
			targS2D2 = 0;
			//NOTE this time gap loops overyear - i.e. Nov->Mar
			//We formulate targS1D2 here from C1 drought plan (publicly available)
			if (doy >= myFex->startDoyMaxR1D2 || doy <= myFex->endDoyMaxR1D2)
			{
				targS1D2 = 1;
			}
			else
			{
				int durationNormRelease = myFex->startDoyMaxR1D2-myFex->endDoyMaxR1D2-1;
				int durationMaxRelease = DAYS_IN_YEAR - durationNormRelease;
				double normRelease = (myFex->linkAnnualCapS1D2 - myFex->linkCapS1D2*durationMaxRelease)/durationNormRelease;
				targS1D2 = normRelease/myFex->linkCapS1D2;
			}
			targR1S1 = targetReleases[0];
			targS1R1 = targetReleases[1];
		}
		free(targetReleases);
		
		//Apply reliability
		if (storageS2 < storageThresholdS2*reliabilityMod)
		{
			targS2D2 = 0;
		}
		if (storageS1 < storageThresholdS1*reliabilityMod)
		{
			targS1D2 = 0;
			targS1R1 = 0;
		}
		/// -------------------------------------------- ///
		
		
		/// ---- CREATE FLOWS FROM TARGET DECISIONS ---- ///
				
		volumeS2D2 = fmax(targS2D2*myFex->linkCapS2D2,myFex->minS2D2);
		if (pipeCondition[0] > 0)
		{
			volumeS2D2 = 0;
		}
		volumeS2D2 = fmin(volumeS2D2,storageS2);

		volumeS1D2 = fmax(targS1D2*myFex->linkCapS1D2,myFex->minS1D2);
		if (pipeCondition[1] > 0)
		{
			volumeS1D2 = 0;
		}

		//NOTE this time gap loops within year - i.e. Apr->Nov
		if (doy >= myFex->startDoyNoR1S1 && doy <= myFex->endDoyNoR1S1)
		{
			volumeR1S1 = 0;
		}
		else
		{
			volumeR1S1 = targR1S1*myFex->linkCapR1S1;
			if (pipeCondition[2] > 0)
			{
				volumeR1S1 = 0;
			}
			
			volumeR1S1 = fmin(volumeR1S1,fmax((fexFlowR1 - myFex->R1Thresh)/2,0)); //where myFex->R1Thresh = (1.16*86.4)Ml/d
			volumeR1S1 = fmin(volumeR1S1,fmax((fexFlowR2 + fexFlowR1 - myFex->R2Thresh),0)); //where myFex->R2Thresh = (3.158*86.4)Ml/d
		}

		double waterAvailableAtR2 = fmax((fexFlowR2 + fexFlowR1 - myFex->R2Thresh),0); //How much water you can take without breaking R2 license (any further abstraction will have to be released from S1)
		double waterAvailableAtR2_ = fexFlowR1*myFex->r2r1ConversionFactor; 
		double waterAvailableForC1 = fmin(waterAvailableAtR2_,waterAvailableAtR2);
		double maximumReleaseWithoutWaste = fmax(myFex->fexDemandD1Fixed - waterAvailableForC1,0);
		
		volumeS1R1 = fmax(targS1R1*myFex->linkCapS1R1,0);
		volumeS1R1 = fmin(volumeS1R1,maximumReleaseWithoutWaste);

		/// -------------------------------------------- ///
		
		
		/// ------------ RESOLVE CONSTRAINTS ----------- ///
		//Check that license allows abstraction
		volumeS1R1 = fmin(volumeS1R1,fmax(myFex->linkAnnualCapR1S1 - S1R1AnnualTotal,0));
		volumeR1S1 = fmin(volumeR1S1,fmax(myFex->linkAnnualCapR1S1 - R1S1AnnualTotal,0));
		//No simulataneous release and abstraction from R1
		if (volumeR1S1 > 0 && volumeS1R1 > 0)
		{
			if (volumeR1S1 > volumeS1R1)
			{
				volumeR1S1 = volumeR1S1 - volumeS1R1;
				volumeS1R1 = 0;
			}
			else
			{
				volumeS1R1 = volumeS1R1 - volumeR1S1;
				volumeR1S1 = 0;
			}
		}

		//Mass balance on s1leball from R1
		storageS1 = storageS1 + volumeR1S1;
		volumeS1R1 = fmin(volumeS1R1,storageS1); // we currently draw s1R1 followed by s1D2 (on the basis that C1 own the dam)
		storageS1 = storageS1 - volumeS1R1; 
		if (myFex->simulationType !=1)
		{
			volumeS1D2 = fmin(volumeS1D2,storageS1);
		}

		//Check that license allows abstraction
		volumeS1D2 = fmin(volumeS1D2,fmax(myFex->linkAnnualCapS1D2 - S1D2AnnualTotal,0)); 
		volumeS2D2 = fmin(volumeS2D2,fmax(myFex->linkAnnualCapS2D2 - S2D2AnnualTotal,0)); 
		volumeS2D2 = fmin(volumeS2D2,fmax(myFex->linkMonthlyCapS2D2 - S2D2MonthlyTotal,0)); 

		//Find supply to D2
		supplyToD2 = volumeS1D2 + volumeS2D2;
		if (supplyToD2 > fexDemandD2 && myFex->simulationType != 2)
		{
			// Ensure no oversupply of D2
			// we currently reduce oversupply by reducing flows in proportion to the flow size
			
			// we also apply this after the minimum flow requirements (e.g. myFex->minS2D2)
			// so the flows may drop below. You can move this up, but then oversupply may happen
			
			
			
			double err = supplyToD2 - fexDemandD2;
			volumeS2D2 = volumeS2D2 - err*volumeS2D2/supplyToD2;
			volumeS1D2 = volumeS1D2 - err*volumeS1D2/supplyToD2;
			supplyToD2 = supplyToD2 - err;
		}

		//Remove d2down supplies from reservoirs
		storageS2 = storageS2 - volumeS2D2;
		storageS1 = storageS1 - volumeS1D2;

		//Update annual counters
		S1D2AnnualTotal = S1D2AnnualTotal + volumeS1D2;
		S2D2AnnualTotal = S2D2AnnualTotal + volumeS2D2;
		S2D2MonthlyTotal = S2D2MonthlyTotal + volumeS2D2;
		R1S1AnnualTotal = R1S1AnnualTotal + volumeR1S1;
		S1R1AnnualTotal = S1R1AnnualTotal + volumeS1R1;
		
		//Calculate C1 supply
		supplyToC1 = fmin(myFex->fexDemandD1Fixed, waterAvailableForC1 + volumeS1R1 - volumeR1S1); 
		volumeBorhole = fmin(myFex->fexDemandD1Fixed - supplyToC1,myFex->linkCapBorehole);
		if (pipeCondition[4] > 0)
		{
			volumeBorhole = 0;
		}		
		supplyToC1 = supplyToC1 + volumeBorhole;


		//Reset annual/monthly counters
		if (doy == DAYS_IN_YEAR)
		{
			R1S1AnnualTotal = 0;
			S1D2AnnualTotal = 0;
			S2D2AnnualTotal = 0;
			S1R1AnnualTotal = 0;
		}
		if ((t+1)%DAYS_IN_MONTH == 0)
		{
			S2D2MonthlyTotal = 0;
		}

		//apply spill
		spillS2 = fmax(0, storageS2 - myFex->storCapS2);
		storageS2 = storageS2 - spillS2;

		spillS1 = fmax(0, storageS1 - myFex->storCapS1);
		storageS1 = storageS1 - spillS1;

		/// -------------------------------------------- ///
		
		
		/// --------- CALCULATE TIME-STEP COSTS -------- ///
		if (t > myFex->warmup)
		{
			gsLinC1 = gsLinC1 + volumeBorhole;
			gsSqrC1 = gsSqrC1 + volumeBorhole*volumeBorhole;
			deficitSqrC1 = deficitSqrC1 + (myFex->fexDemandD1Fixed - supplyToC1)*(myFex->fexDemandD1Fixed    - supplyToC1);
			deficitLinC1 = deficitLinC1 + (myFex->fexDemandD1Fixed  - supplyToC1);
			deficitSqrD2 = deficitSqrD2 + (fexDemandD2 - supplyToD2)*(fexDemandD2 - supplyToD2);
			deficitLinD2 = deficitLinD2 + (fexDemandD2 - supplyToD2);
			if (storageS2 < myFex->storageReliabilityThresholdS2)
			{
				storageReliabilityS2++;
			}
			if (storageS1 < myFex->storageReliabilityThresholdS1)
			{
				storageReliabilityS1++;
			}
			
			NPVc2 = NPVc2 + (volumeS1D2*myFex->linkCostS1D2)/pow((1+myFex->discountRate),t);
			NPVC1 = NPVC1 + (volumeR1S1*myFex->linkCostR1S1)/pow((1+myFex->discountRate),t);
			NPVtotal = NPVc2 + NPVC1;
			pumpCostC2 = pumpCostC2 + volumeS1D2*myFex->linkCostS1D2;
			pumpCostC1 = pumpCostC1 + volumeR1S1*myFex->linkCostR1S1 + volumeBorhole*myFex->linkCostGroundBorehole;
		}
		/// -------------------------------------------- ///
		
		
		/// -------------- PRINT DYNAMICS -------------- ///
		if (myFex->printDynamics == 1)
		{
			if (t == 0)
			{
				fprintf(fid,"%f %f ", myFex->fexInitialStorS2,myFex->fexInitialStorS1);
			}
			else
			{
				fprintf(fid,"%f %f ", storageS2, storageS1);
			}
			fprintf(fid,"%f %f %f %f ", fexFlowS2, fexFlowS1, fexFlowR1,fexFlowR2);
			fprintf(fid,"%f %f ",fexDemandD2,myFex->fexDemandD1Fixed);
			fprintf(fid,"%f %f %f %f ",volumeS2D2,volumeS1D2,volumeR1S1,volumeS1R1);
			fprintf(fid,"%f %f %f %d ",volumeBorhole,spillS2,spillS1,doy);
			fprintf(fid,"%f %f %f %f\n",evapS2, evapS1,s2Env, s1Env);
			
		}
		
		/// -------------------------------------------- ///
		if (storageS2 < 0)
		{
			if ((storageS2 < -(DOUBLE_PRECISION*10000)) && myFex->simulationType != 2)
			{
				printf("warning: s2 storage below 0");
			}
			storageS2 = 0;
		}
		
		if (storageS1 < 0)
		{
			if ((storageS1 < -(DOUBLE_PRECISION*10000)) && myFex->simulationType != 1)
			{
				printf("warning: s1 storage below 0");
			}
			storageS1 = 0;
		}
	}
	
	//edit: LoRam 4/Apr/2018
	for (i=0; i < myFex->parent->NUM_AUTOCORRELATED_PROCESSES;i ++)
	{
		free(innovationHolders[i]);
		free(residualHolders[i]);
	}
	free(innovationHolders);
	free(residualHolders);	
	free(pipeCondition);
	free(timeUntilBreak);
	free(breakLambdas);
	free(allParma);

	if (myFex->simulationType == 3)
	{
		free(c2Weights);
		free(C1Weights);
	}
	
	
	/// ----------- CALCULATE OBJECTIVES ----------- ///
	
	//THIS MAKES ANGELS WEEP... YOU NEED AN OBJECTIVE SPECIFIER INDEX
	
	l = 0;
	if (myFex->objectiveSpecifier[0] == 1)
	{
		// timeObjectives[l] = storageReliabilityS2;
		objs[l] = storageReliabilityS2/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(storageReliabilityS2);
	// }
	
	if (myFex->objectiveSpecifier[1] == 1)
	{
		// timeObjectives[l] = storageReliabilityS1;
		objs[l] = storageReliabilityS1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(storageReliabilityS1);
	// }
	
	if (myFex->objectiveSpecifier[2] == 1)
	{
		// timeObjectives[l] = deficitLinD2;
		objs[l] = deficitLinD2/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitLinD2);
	// }
	
	if (myFex->objectiveSpecifier[3] == 1)
	{
		// timeObjectives[l] = deficitSqrD2;
		objs[l] = deficitSqrD2/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitSqrD2);
	// }
	
	if (myFex->objectiveSpecifier[4] == 1)
	{
		// timeObjectives[l] = deficitLinC1;
		objs[l] = deficitLinC1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitLinC1);
	// }
	
	if (myFex->objectiveSpecifier[5] == 1)
	{
		// timeObjectives[l] = deficitSqrC1;
		objs[l] =deficitSqrC1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitSqrC1);
	// }
	
	if (myFex->objectiveSpecifier[6] == 1)
	{
		// timeObjectives[l] = NPVtotal;
		objs[l] = NPVc2 + NPVC1;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVtotal);
	// }
	
	if (myFex->objectiveSpecifier[7] == 1)
	{
		// timeObjectives[l] = NPVc2;
		objs[l] = NPVc2;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVc2);
	// }
	
	if (myFex->objectiveSpecifier[8] == 1)
	{
		// timeObjectives[l] = NPVC1;
		objs[l] = NPVC1;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVC1);
	// }
	
	if (myFex->objectiveSpecifier[9] == 1)
	{
		// timeObjectives[l] = gsLinC1;
		objs[l] = gsLinC1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(gsLinC1);
	// }
	
	if (myFex->objectiveSpecifier[10] == 1)
	{
		// timeObjectives[l] = gsSqrC1;
		objs[l] = gsSqrC1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(gsSqrC1);
	// }
	
	if (myFex->objectiveSpecifier[11] == 1)
	{
		// timeObjectives[l] = pumpCostTotal;
		objs[l] = pumpCostTotal/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostTotal);
	// }
	
	if (myFex->objectiveSpecifier[12] == 1)
	{
		// timeObjectives[l] = pumpCostC2;
		objs[l] = pumpCostC2/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostC1);
	// }
	
	if (myFex->objectiveSpecifier[13] == 1)
	{
		// timeObjectives[l] = pumpCostC1;
		objs[l] = pumpCostC1/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostC1);
	// }
	if (myFex->printDynamics == 1)
	{
		fprintf(fid,"#\n//Final Storages, S2 S1:\n");
		fprintf(fid,"%f %f\n", storageS2, storageS1);
		
		fprintf(fid,"#\n//Objectives:\n");
		for (i=0;i<myFex->numObjectives;i++)
		{
			fprintf(fid,"%f ",objs[i]);
		}
		fclose(fid);
	}

	
	
	/// -------------------------------------------- ///
}

void sortDescending(int N, double *vec)
{
	//Sort the vector (vec) of doubles (length = N) with largest first -- note that sorted values are stored in the original vector
	int i,j;
	double a;
	for (i = 0; i < N; ++i)
    {
        for (j = i + 1; j < N; ++j)
        {
            if (vec[i] < vec[j])
            {
                a =  vec[i];
                vec[i] = vec[j];
                vec[j] = a;
            }
        }
    }
}

void sortAscending(int N, double *vec)
{
	//Sort the vector (vec) of doubles (length = N) with smallest first -- note that sorted values are stored in the original vector
	int i,j;
	double a;
	for (i = 0; i < N; ++i)
    {
        for (j = i + 1; j < N; ++j)
        {
            if (vec[i] > vec[j])
            {
                a =  vec[i];
                vec[i] = vec[j];
                vec[j] = a;
            }
        }
    }
}

void calcPercentile(double *percentiles, double *ensembleObjectives, struct fexEnsemble *myEnsemble,int numObjectives)
{
	//Implementation of MATLAB prctile function
	double pctile = (myEnsemble->percentile/100)*myEnsemble->ensembleSize;
	int index = (int)(pctile+0.5);
	double remainder = pctile - index;
	int i, j;
	for (j=0;j< numObjectives;j++)
	{
		double *ensembleSingleObjective = calloc(myEnsemble->ensembleSize,sizeof(double));
		for (i=0;i<myEnsemble->ensembleSize;i++)
		{
			ensembleSingleObjective[i] = ensembleObjectives[j + i*numObjectives];
		}
		sortAscending(myEnsemble->ensembleSize,ensembleSingleObjective);
		if (index < 1)
		{
			percentiles[j] = ensembleSingleObjective[0];
		}
		else if (index >= myEnsemble->ensembleSize)
		{
			percentiles[j] = ensembleSingleObjective[myEnsemble->ensembleSize - 1];
		}
		else
		{
			percentiles[j] = (0.5 - remainder)*ensembleSingleObjective[index-1] + (0.5 + remainder)*ensembleSingleObjective[index];
		}
		free(ensembleSingleObjective);
	}
}

void calcMean(double *means, double *ensembleObjectives, struct fexEnsemble *myEnsemble,int numObjectives)
{
	double runsum;
	//Implementation of MATLAB prctile function
	int i, j;
	for (j=0;j< numObjectives;j++)
	{
		runsum = 0;
		for (i=0;i<myEnsemble->ensembleSize;i++)
		{
			runsum = runsum + ensembleObjectives[j + i*numObjectives];
		}
		means[j] = runsum/myEnsemble->ensembleSize;
	}
}
void evaluateEnsembleMean(double* weights, double* objs, double* consts, struct fexEnsemble *myEnsemble)
{

	//A wrapper around simFex that evaluates a fexEnsemble and returns the percentiles in objs.
	int i,j,numObjectives = myEnsemble->ensemble[0]->numObjectives;
	double *ensembleObjectives = calloc(myEnsemble->ensembleSize*numObjectives,sizeof(double));
	FILE *f;
	if (myEnsemble->ensemble[0]->printDynamics == 1)
	{
		 f = fopen(myEnsemble->ensemble[0]->dynamicFid,"w");
	}
	for (i=0;i< myEnsemble->ensembleSize;i++)
	{
		double *holderObjectives = calloc(numObjectives,sizeof(double));
		
		simFex(weights, holderObjectives, 0,myEnsemble->ensemble[i]);
	
		for (j = 0;j < numObjectives; j++)
		{
			if (myEnsemble->ensemble[i]->printDynamics == 1)
			{
				fprintf(f,"%f",holderObjectives[j]);
				if (j == numObjectives -1)
				{
					fprintf(f,"\n");
				}
				else
				{
					fprintf(f," ");
				}
			}
			ensembleObjectives[j + i*numObjectives]= holderObjectives[j];
		}
		free(holderObjectives);
	}
	if (myEnsemble->ensemble[0]->printDynamics == 1)
	{
		fprintf(f,"#\n");
		fclose(f);
	}
	// calcPercentile(objs, ensembleObjectives, myEnsemble, numObjectives);
	double *means = calloc(numObjectives,sizeof(double));
	calcMean(means, ensembleObjectives, myEnsemble, numObjectives);
	for (j = 0; j < numObjectives;j++)
	{
		objs[j] = means[j];
	}
	free(ensembleObjectives);
}

struct fexEnsemble createMultiFramingEnsemble(int numFramings, int* framingsIncluded, int ensembleSizePerFraming, int ensembleT, int baseSeed, double percentile, int simType)
{
	int i, j;
	struct fexEnsemble combinedEnsemble;
	combinedEnsemble.ensembleSize = ensembleSizePerFraming*numFramings;
	combinedEnsemble.ensembleT = ensembleT;
	combinedEnsemble.baseSeed = baseSeed;
	combinedEnsemble.percentile = percentile;
	combinedEnsemble.simulationType = simType;		
	createEnsemble(&combinedEnsemble,framingsIncluded[0]);
	srand(baseSeed);
	for (i = 0; i < numFramings; i++)
	{
		struct fexEnsemble *myEnsemble = malloc(sizeof(struct fexEnsemble));
		myEnsemble->ensembleSize = 0;
		myEnsemble->ensembleT = 0;
		myEnsemble->baseSeed = rand();
		myEnsemble->percentile = percentile;
		myEnsemble->simulationType = simType;
		createEnsemble(myEnsemble,framingsIncluded[i]);	
		for (j = 0; j < ensembleSizePerFraming; j++)
		{
			combinedEnsemble.ensemble[i*ensembleSizePerFraming + j]->parent = myEnsemble;
		}
		// destroyEnsemble(&myEnsemble);
	}
	return combinedEnsemble;
}

/*
void repopulateFex(struct fexProblem *myFex, struct fexEnsemble *myEnsemble,int seed)
{
	destroyFex(myFex);
	createFex(myFex, myEnsemble, myEnsemble->ensembleT, seed, myEnsemble->simulationType);
	populateSeededVariables(myFex, myEnsemble, seed);
}*/
///	---------------- EXAMPLE SIMULATION ---------------- ///
// compiles (on cygwin) with command gcc -o fex.exe fex.c qr_solve.c r8lib.c -lm

// int main(){
	// struct fexEnsemble myEnsemble;
	// myEnsemble.ensembleSize = 100;
	// myEnsemble.ensembleT = 1000;
	// myEnsemble.baseSeed = 1;
	// myEnsemble.percentile = 99;
	// myEnsemble.simulationType = 0;
	
	// createEnsemble(&myEnsemble,6);
		
	// int i;
	// double *weights = calloc(myEnsemble.ensemble[0].cooperativeNet->numWeights,sizeof(double));
	// for (i = 0; i < myEnsemble.ensemble[0].cooperativeNet->numWeights;i++)
	// {
		// weights[i] = 0;
	// }
	// double *objective = calloc(myEnsemble.ensemble[0].numObjectives,sizeof(double));

	// evaluateEnsemble(weights, objective, 0, &myEnsemble);

	// for (i = 0; i < myEnsemble.ensemble[0].numObjectives;i++)
	// {
		// printf("%f ",objective[i]);
	// }
	// printf("\n");
	// destroyEnsemble(&myEnsemble);
	// free(objective);
	// free(weights);
// }

