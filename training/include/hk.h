#ifndef HK_H
#define HK_H
#include <stsp.h>
#include <node.h>
#include <heap.h>
#include <vector>
#include <limits>
#include <exception>
#include <stdexcept>
#include <string>



class parentMaxEdge {

	public :
	int parentIndex;	/*!< Ancestor index (an ancestor is a node with a degre >2 in the spanning tree).*/
	int maxEdgeSourceIndex;	/*!< Index of the source of the maximum edge in the spanning tree on the path to ancestor.*/
	int branche;	/*!< Index of ancestor's subtree this node is into.*/
	double weight;	/*!< Weight the maximum edge in the spanning tree on the path to ancestor.*/
	
	parentMaxEdge() : 
	parentIndex(-1),
	maxEdgeSourceIndex(-1),
	weight(-1.0),
	branche(-1)
	{
	};
};
class HeldKarp {

	SymTSP* graphe;  /*!< Pointer to the TSP graphe.*/
	//SymTSP* min1tree;
	
	//int MaxEdgeMieux;
	//int MaxEdge;
	
	//char * euh;
	
	int dfscount;
	

	Pointeur * oneNode;		/*!< Pointer to the 1-node of the 1-tree.*/
	
	Pointeur * mstRoot;		/*!< Pointer to the root of the directed minimum spanning tree computed with Prim.*/
	




	std::vector<double> penalite; /*!< Penalties on nodes.*/
	std::vector<double> bestPenalites; 
	
	std::vector<int> degrees; 
	std::vector<int> bestDegrees; 

	std::vector<bool> nodeExist; /*!< True if node i exist in graphe.*/


//	Pointeur* pointeurs[graphe->getMaxSize()];

	//std::vector<Edge*> oneTreeEdges;


	std::vector<parentMaxEdge> parents;		 /*!< Ancestors of each node.*/
	
	double oneTreeWeight;		/*!< Weight of the 1-tree.*/
	
	double HKBound;				/*!< Held and Karp bound value.*/
	bool computed;				/*!< True if the Held and Karp bound has been computed.*/
	
	
	/*!
	*  \brief Return penalty value on the edge from source to dest.
	*/
	double getPenalite(int source, int dest);
	
	/*!
	*  \brief Return the weight of the heaviest edge on the cycle created when adding edge from source to dest
	* to the spanning tree of the 1-tree.<br>
	*
	*This method find the cycle created by adding edge (source,dest) to the minimum spanning tree 
	*using the directed tree structure given by Prim algorithm.<br>
	*The following pseudo-code is used :<br>
	*\verbatim
	*p[v] : father of v in the directed tree
	*d[v] : weight of the edge connecting v to p[v]
	*rank[v] : rank of node v in the spanning tree
	*maxWeight=-infinty
	*
	*i=source
	*j=dest
	*
	*if( rank[j] > rank[i] ) {
	*	swap(i,j)
	*}
	*
	*while(rank[i] > rank[j]) {
	*	maxWeight = max( d[i], maxWeight)
	*	i=p[i]
	*}
	*
	*while(i != j) {
	*	maxWeight = max( d[i], maxWeight)
	*	i=p[i]
	*	maxWeight = max( d[j], maxWeight)
	*	j=p[j]
	*}
	*return maxWeight
	*\endverbatim
	*<br><br>
	*Expected complexity O(n)
	*/
	void getWeightOfMaxEdgeOnCycleCreatedBy(Edge* edge, bool * edgeIsInOneTree, double * weightOfRemovedEdge);
	
	
	
	/*!
	*  \brief  Return the weight of the heaviest edge on the cycle created when adding edge from source to dest
	* to the spanning tree of the 1-tree. Use the directed spanning tree to speed up computation.<br>
	*
	* After a call to dfs(), each node knows its ancestor and the heaviest edge on the path to its ancestor.<br>
	* The cycle created by adding edge (source,dest) to the spanning tree can be efficiently found by comparing 
	* ancestors of source and dest and jumping from ancestor to ancestor.<br>
	*
	* This pseudo-code is used :
	* \verbatim
	*p[v] : father of v in the directed tree
	*d[v] : weight of the edge connecting v to p[v]
	*a[v] : ancestor of v
	*w[v] : weight of heaviest edge on path from v to its ancestor
	*rank[v] : rank of node v in the spanning tree
	*maxWeight=-infinty
	*
	*
	*i=source
	*
	*#We move up until both nodes has the same ancestor
	*while( a[i] != a[j] ) {
	*	if(rank[i] > rank[j]) {
	*		maxWeight= max( w[i], maxWeight)
	*		i=a[i]
	*	}
	*	else
	*		maxWeight= max( w[j], maxWeight)
	*		j=a[j]
	*	}
	*}
	*
	*if(i==j)
	*	return maxWeight
	*
	*if( rank[j] > rank[i] ) {
	*	swap(i,j)
	*}
	*
	*if(i and j are not in the same subtree) {
	*	maxWeight = max( w[i], w[j], maxWeight);
	*	return maxWeight;
	*}
	*
	*if(i and j are in the same subtree)  {#then j is in the path from i to a[i]
	*	while(i != j) {
	*		maxWeight= max( d[i], maxWeight)
	*		i=p[i]
	*	}
	*	return maxWeight;
	*}
	*\endverbatim
	*<br><br>
	*/
	pair<double,bool> getWeightOfMaxEdgeOnCycleCreatedBy2(Edge* ed);
	
	

	
	
	/*!
	*  \brief Return the weight of the edge to remove from the 1-tree if the edge (source,dest) is forced.<br>
	*
	*Check if edge (source,dest) is connected to the 1-node. If it is not, call getWeightOfMaxEdgeOnCycleCreatedBy(int source, int dest).
	*/
	void getWeightOfRemovedEdge(Edge* edge, bool * edgeIsInOneTree, double * weightOfRemovedEdge);
	

	
	

	
	
	/*!
	*  \brief Compute ancestors of each node by doing a depth first search in the minimum spanning tree.<br>
	*
	* An ancestor of node x is the first node with degre >2 found when going from x toward to thee in the directed spanning tree
	* obtained with Prim algorithm.<br>
	*
	* Ancestors for each node and the heaviest edge on the path to ancestor can easily be found by performing a depth first search from the root of the tree. <br>
	*
	* This pseudo-code if used :<br>
	*\verbatim
	*p[v] : father of v in the directed tree
	*d[v] : weight of the edge connecting v to p[v]
	*a[v] : ancestor of v
	*w[v] : weight of heaviest edge on path from v to its ancestor

	*dfs(x, a,w) {
	*	a[x]=a
	*	if(d[v] > w) {
	*		w=d[v]
	*	}
	*	w[v]=w
	*
	*	if(x is a leaf) STOP
	*	if( degre[x] > 2 ) {
	*		a=x
	*		w=-infinity
	*	}
	*	for( y succesor of x in the spanning tree) {
	*		dfs(x,a,e)
	*	}
	*
	*}
	*\endverbatim
	*<br><br>
	*Expected complexity O(n)
	*/
	void dfs(Pointeur* source, parentMaxEdge pme);
	
	
	
	public :
	HeldKarp(SymTSP * g);
	
	
	std::vector<Pointeur*> pointeurs; /*!< Pointeurs used in Prim. Also contains the MST*/
	
	double oneEdge1Weight, oneEdge2Weight, oneEdge3Weight; /*!< Weights of edges connected to the 1-node in the 1-tree.*/
	
	
	int index1node;			/*!< Index of the 1-node of the 1-tree.*/
	
	Edge * oneEdge1;
	Edge * oneEdge2;	/*!< Edges connected to the 1-node in the 1-tree.*/
	int index1Edge1, index1Edge2;
	/*!
	*  \brief Compute Held and Karp bound.<br>
	*
	* Increments penalties using the following scheme :
	*\verbatim
	*alpha=2
	*
	* for(int l = 1; l < 25; l++ ) {
	*	for(int k=1; k< 0.15*n; k++) {
	*		computeMinimumOneTree
	*		HKbound = oneTreeWeight - 2*sumOfPenalties
	*		violation = sum( degre[i] -2 )
	*		step = alpha*(upperBound - HKBound)/violation
	*		if(step < 0.01)
	*			STOP
	*		penalty[i] +=step*(degre[i] - 2)
	*
	*	}
	*	alpha=0.5*alpha
	*}
	*\endverbatim
	*/
	void computeBound(double upperBound, bool *tourFound, bool *canStopBranching, float *HKbound_fnf, float *factors);
	
		/*!
	*  \brief Compute the minimum 1-tree.<br>
	*
	* Use Prim algorithm, using custom made Heap, on the graphe where the 1-node has been removed.<br>
	* The two minimum edges connected to the 1-node are the added.
	*/
	void PrimMieux(bool * canContainTour);
	
	/*!
	*  \brief Compute Held and Karp bound with a different scheme.<br>
	*
	* The only difference with computeBound() is that when the bound lower, 
	* this method reduce alpha (thus the step size).<br><br>
	*
	* The convergence of this scheme is faster thant computeBound(), but it leads to weaker bounds.<br>
	*
	* Increments penalties using the following scheme :<br>
	*\verbatim
	*alpha=2
	*computeMinimumOneTree()
	*HKbound = oneTreeWeight - 2*sumOfPenalties
	* for(int l = 1; l < 500; l++ ) {
	*	while(HKBound > oldBound) {
	*		oldBound=HKBound
	*		computeMinimumOneTree()
	*		HKbound = oneTreeWeight - 2*sumOfPenalties
	*		violation = sum( degre[i] -2 )
	*		step = alpha*(upperBound - HKBound)/violation
	*		if(step < 0.01)
	*			STOP
	*		penalty[i] +=step*(degre[i] - 2)
	*
	*	}
	*	alpha=0.5*alpha
	*}
	*\endverbatim
	*/

	void filter_reducedCosts(double upperBound,std::vector<Edge*>* toForce);
	
	void filter_forcedDegrees( std::vector<Edge*>* toRemove);
	
	
	std::vector<Edge*>* filter3(double upperBound);
	pair<double,bool> getWeightOfMaxEdgeOnCycleCreatedBy3(int source, int dest, double edgeW);
	pair<double,bool> getWeightOfRemovedEdge3(int source, int dest, double edgeW);

	void clearPenalties();
	
	double getBound();
	
	bool isEdgeInMST(int source, int dest);
	void printDot(string nomFichier);
		void printDot2(string nomFichier);
	void printGrapheWithPenalties(string nomFichier);
	
	
	std::vector<double> getPenalites();
	
	void setPenalites(std::vector<double> &vecteur);
	
	double getOneTreeWeightSansPenalite();
	
	bool isInOneTree(int source, int dest) const;
	
	void PrimForce(bool * containsTour, std::vector<Edge*>* toForce);
	
	void isTour(bool * isTour, bool * canContainTour);
	
	void forceMSTEdges(double upperBound, std::vector<Edge*>* toForce);
	
	void forceEdgesWithDegree2(std::vector<Edge*>* toForce);
	
	std::vector<Edge*> toRemove;
	
	bool filterOK;
	
	double getSumOfCostOfMinEdgesOnLeaves();
	
	bool isNodeWithForcedDegree3();
	
	bool isNodeWithDegree1();
	pair<double,bool> getWeightOfRemovedEdgePrint(Edge* ed, bool mieux, ofstream fb);
	
	void filterPrint();
	pair<double,bool> getWeightOfMaxEdgeOnCycleCreatedByPrint(Edge* ed, ofstream fb);

};
#endif
