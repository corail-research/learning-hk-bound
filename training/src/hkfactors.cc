#include <hkfactors.h>


/* global parameters */
const char* FileIn;
const char* FileOut;
const char* FileSol;
const char* dotFile;
bool isSym;
bool tour_file_specified;
bool write_output;
int output_level;
bool prim;
bool prim2;
bool kruskal;
bool printDot;
double upperBound;
bool upperBoundSpecified;
bool heldKarp;
bool filter;
bool heldKarp2;
bool filter2;
bool forceCost;
bool forceCut;
int nbForce;
bool heavyFiltering;
bool noSecondFiltering;
bool noFiltering;
bool forceDegree;
bool filterDegree;

bool printReducedGraph;
const char* reduceGraphFileName;

bool useAPBound;
bool filterAP;

bool ap_bb;

bool useReducedCostInAP;

bool verbose = false;



extern "C"{

float *hkfactors(int argc, char* argv[]){

cout << "Start HK bound" << endl;
/*cout << argc << endl;
cout << argv[0] << endl;
cout << argv[1] << endl;
cout << argv[2] << endl;
cout << argv[3] << endl;
cout << argv[4] << endl;
cout << argv[5] << endl;
cout << argv[6] << endl;
cout << argv[7] << endl;
cout << argv[8] << endl;
cout << argv[9] << endl;
cout << argv[10] << endl;
cout << argv[11] << endl;
cout << argv[12] << endl;
cout << "printed" << endl;*/
	//return argc;}

	parseArgs(argc, argv);
	
	//upperBound+=0.1;
	
	int n=3;				// number of cities

	//readTSPInstance(n, YOUR_DATA_STRUCTURE);
	//Graphe* g =readTSPInstance(n);

	Graphe* g = readInputFile();
	//Graphe* g = readInputMatrix();
	if(verbose) cout <<"Parsing done."<<endl;
	
	if(printDot) {
		g->printDot(dotFile);
	}
	SymTSP* stsp;
	if(!isSym) {
	    	
		ATSP* atsp =dynamic_cast<ATSP*>(g);
		stsp = atsp->getTSP();	
	}
	
	else {
		stsp =dynamic_cast<SymTSP*>(g);
	}	

	cout <<"Number of nodes : "<<stsp->getSize()<<endl;
	cout << "Number of edges : "<< (*stsp->getEdges()).size()<<endl;

	if( !upperBoundSpecified) {
		cout<<"Need an upper bound to launch branch and bound."<<endl
		<<"Please specify one using -UB ."<<endl;
		exit(0);
	}
	

	int size_bf = stsp->getSize()+1;
	float *bound_factors = new float[size_bf];
	//float *bound_factors = malloc(sizeof(float)*size_bf);
	
	
	BB bb = BB(stsp, upperBound, nbForce);
	if(ap_bb) {
		cerr<<"Lauching AP based branch and bound."<<endl;
		useReducedCostInAP=false;
		bb.compute_ap();
	}
	else {
		cerr<<"Lauching HK based branch and bound."<<endl;
		useReducedCostInAP=true;
		//bb.compute();
		bb.compute(size_bf, bound_factors); //bound_factors = [HKbound,penality_factors]
	}

	cerr << "HK BOUND" << endl;
 	cerr << bound_factors[0] << endl;
	for(int i=1;i<size_bf;i++){
		cerr << bound_factors[i] << " ";
	}
	
	
	//BiGraphe bg(stsp);
	
	//bg.print(); 


  return bound_factors;
}

void hkfact_del(float *bound_factors){delete [] bound_factors;}

}
