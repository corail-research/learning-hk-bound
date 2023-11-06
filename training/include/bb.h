#ifndef BB_H
#define BB_H

#include <hk.h>
#include <stsp.h>
#include <hungarianMethod.h>

#include <time.h>

class BB {

	private :
	HeldKarp hk;
	SymTSP * graphe;
	
	hungarianMethod hg;
	
	clock_t startTime;
	
	float percent_edges_filtered;
	
	double upperBound;
	int nbNode;
	bool bb_tourFound;
	double tourCost;
	
	
	std::vector<Edge*>* toForce_temp;
	std::vector<Edge*>* toRemove_temp;
	
	pair<bool,Edge*> getEdgeWithMaxReplacementCost();
	void getEdgesToBranchOn(std::vector<Edge*>* edgesToBranchOn);
	
	int nbForce;
	
	int nbForceCut;
	int nbForceCost;
	
	public :
	BB(SymTSP * graphe, double UB, int nF);
	
	~BB();

	void test();

	void compute(int size_bf, float bound_factors[]);

	
	void dfs_Remove(int profondeur);
	
	void filter_and_force(std::vector<Edge*>* toRemove, std::vector<Edge*>* toForce, bool *canStopBranching, float *HKbound_fnf, float *factors);
	
	
	void compute_ap();

	
	void dfs_Remove_ap(int profondeur);
	
	void filter_and_force_ap(std::vector<Edge*>* toRemove, std::vector<Edge*>* toForce, bool *canStopBranching);
	void printEndInfos(bool isOptimal);

};

#endif 
