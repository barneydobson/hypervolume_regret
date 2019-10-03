#ifndef FEX_H_INCLUDED
#define FEX_H_INCLUDED

struct fexProblem {
	//Hard constraints
	double linkCapS2D2;
	double linkCapS1D2;
	double linkCapR1S1;
	double linkCapS1R1;
	double linkCapBorehole;
	
	double linkAnnualCapS1D2;
	double linkAnnualCapR1S1;
	double linkAnnualCapS2D2;
	double linkMonthlyCapS2D2;
	double linkAnnualCapS1R1;
	
	double storCapS2;
	double storCapS1;
	
	double discountRate;	
	
	double deadStorS2;
	double deadStorS1;
	
	int startDoyNoR1S1;
	int endDoyNoR1S1;
	
	int startDoyMaxR1D2;
	int endDoyMaxR1D2;
	
	double storElevationS2[2];
	double elevationAreaS2[2];
	
	double storElevationS1[2];
	double elevationAreaS1[2];
	
	double s1CompensationFixed;
	double s2CompensationFixed;
	
	double R1Thresh;
	double R2Thresh;
	
	//Synthetic Generation
	
	
	
	//Forcing
	int fexT;
	
	double *fexFlowsS2;
	double *fexFlowsS1;
	double *fexFlowsR1;
	double *fexFlowsR2;
	
	double r2r1ConversionFactor;
	
	double *s1FishReleases;
	
	int *fexDoys;
	
	double *linkCapS2D2TimeVarying;
	double *linkCapS1D2TimeVarying;
	double *linkCapR1S1TimeVarying;
	double *linkCapS1R1TimeVarying;
	double *linkCapBoreholeTimeVarying;
	
	double *fexPets;
	//double *W_temps;
	
	double *fexDemandsD2;
	// double *fexDemandsSww;
	double fexDemandD1Fixed;
	
	double fexInitialStorS2;
	double fexInitialStorS1;
	
	//Evaluation parameters
	double linkCostS2D2;
	double linkCostS1D2;
	double linkCostR1S1;
	double linkCostS1R1;
	double linkCostGroundBorehole;
	
	double storageReliabilityThresholdS2;
	double storageReliabilityThresholdS1;
	
	int numObjectives;
	int numObjectiveOptions;
	int *objectiveSpecifier;
	// int *objectiveSpecifierIndex;
	
	//Decision model parameters
	double inflowNormMaxMinS2[2];
	double inflowNormMaxMinS1[2];
	double inflowNormMaxMinR1[2];
	double inflowNormMaxMinR2[2];
	
	double minS2D2;
	double minS1D2;
	
	double petNormMaxMinS1[2];
	
	double demandNormMaxMinD2[2];
	
	
	int numInputOptions;
	struct netOptions *cooperativeNet;
	struct netOptions *c2Net;
	struct netOptions *c1Net;
	
	//edit:LoRam 4/Apr/18
	struct fexEnsemble *parent;
	int seed;
	int warmup;
	//endedit
	
	//Simulation settings
	int simulationType;
	int printDynamics;
	char *dynamicFid;
	// int includetimeObjectives;
	// double **timeObjectives;
};

struct netOptions{
	int numNetInputs;
	int *policyInputSpecifier;
	int *policyInputSpecifierIndex;
	
	int numNetOutputs;
	int betaInd;
	int betaIndMax;
	int numLayers;
	int *architecture;
	int numWeights;
};

struct parma{
	int numHarmonics;
	double *harmonics;
	double *thetas;
	
	int numAR;
	double *AR;
	
	int numMA;
	double *MA;
	
	double cons;
	double varc;
};
struct fexEnsemble{
	int ensembleT;
	int ensembleSize;
	int baseSeed;
	int NUM_AUTOCORRELATED_PROCESSES;
	double percentile;
	
	struct fexProblem **ensemble;
	
	//Synthetic Generation parameters
	struct parma S2Parma;
	struct parma S1Parma;
	struct parma R1Parma;
	struct parma R2Parma;
	struct parma D2Parma;
	struct parma petParma;

	struct parma s2Theta;
	struct parma s1Theta;
	
	int offset;
	
	double **cholesky;
	
	double breakLambdaProbabilityS2D2[2];
	double breakLambdaProbabilityS1D2[2];
	double breakLambdaProbabilityR1S1[2];
	double breakLambdaProbabilityS1R1[2];
	double breakLambdaProbabilityBorehole[2];
	
	double s1FisheriesAnnual;
	double maxDailyFishRelease;
	int startAllowableFishReleaseDoy;
	int endAllowableFishReleaseDoy;
	int maxFishReleaseLength;
	int minFishReleaseLength;
	
	int simulationType;
	int framing;
};
void destroyEnsemble(struct fexEnsemble *myEnsemble);
void createEnsemble(struct fexEnsemble *myEnsemble, int framing);

void evaluateEnsembleMean(double* weights, double* objs, double* consts, struct fexEnsemble *myEnsemble);
struct fexEnsemble createMultiFramingEnsemble(int numFramings, int* framingsIncluded, int ensembleSizePerFraming, int ensembleT, int baseSeed, double percentile, int simType);
void evaluateNet(double *outputs, double *inputs, int *architecture, int numLayers, double *weights);
#endif