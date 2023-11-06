


//////////////////////////////////////////////////
#define DBG(x)       // use this to not print debug info
//#define DBG(x) x     // use this to print debug info


#include <definitions.h>
#include <parse.h>
#include <hk.h>
#include <bb.h>
#include <bigraphe.h>
#include <hungarianMethod.h>
//#include <sys/time.h>


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

int main(int argc, char* argv[]){

	parseArgs(argc, argv);
	
	//upperBound+=0.1;
	
	int n=3;				// number of cities

	//readTSPInstance(n, YOUR_DATA_STRUCTURE);
	//Graphe* g =readTSPInstance(n);
	Graphe* g =readInputFile();
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
	float bound_factors[size_bf];
	
	
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

	/*
	Tour* tour;
	if(tour_file_specified) {
		tour = readTourInstance();
		//tour->print();
		double cost = g->getTourValue(*tour);
		cout << "tour cost : "<<cost << endl;
	}
*/
	//g->print();


  return 0;
}
