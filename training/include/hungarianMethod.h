#ifndef HUNG_H
#define HUNG_H
#include <bigraphe.h>
#include <vector>
#include <limits>
#include <hk.h>


class hungarianMethod {

const static bool debug=false;

private:

BiGraphe* bigraphe;
bool isFromATSP;

list<BiNode*> leaves;
list<BiNode*> leftNodes_toMatch;
list<BiNode*> leftNodes_inTree;
list<BiNode*> rightNodes_notInTree;
list<BiNode*> rightNodes_inTree;


int matchingSize;
double matchingCost;




void initPotential(bool * canContainTour);

void init(binodeEdges * toRemove);

BiNode *  findAlternatePath(BiNode * startNode, bool * canContainTour );

BiNode *  expandAlternateTree(BiNode * node, bool * matchingFound);


void switchPath(BiNode * startNode, BiNode * node_endOfPath);


public :

hungarianMethod(SymTSP* stsp);
~hungarianMethod();

void computeBound(bool * canContainTour);

double getBound();


void filter(double upperBound, double lowerBound, std::vector<Edge*>* toRemove);

void printDot(string nomFichier);

void removeEdges(std::vector<Edge*>* toRemove);
void removeEdge(Edge* edge);
void addEdges(std::vector<Edge*>* toRemove);
void addEdge(Edge* edge);

void setEdgesCostsInHK();
void getEdgesToBranchOn(std::vector<Edge*>* edgesToBranchOn);

bool isTour();
};
#endif
