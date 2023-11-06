

#ifndef GRAPHE_H
#define GRAPHE_H


#include <set>
#include <list>
#include <map>
#include <vector>
#include <node.h>
#include <tour.h>
#include <iostream>
#include <string>
#include <fstream>
//#include <stsp.h>
/*!
 * \file Graphe.h
 * \brief Header for Graphe.
 * \author Pascal Benchimol
 */
typedef std::vector<Node*> graphNodes;


using namespace std;

/*!
 * \class Graphe
 * \brief A superclass for SymTSP and ATSP. Represents a graph.
 */
class Graphe {
protected:


int maxSize;
graphNodes nodes; /*!< Set of nodes*/
graphEdges edges;	/*!< Set of edges*/

public:

/*!
 *  \brief Constructeur
 *  \param size : number of nodes.
 */
Graphe(int size);
Graphe();
~Graphe();


/*!
 * \brief Return the maximum index of a node that can be found in this graphe.
 */
int getMaxSize() const;


void removeEdge(Edge * e, bool * canContainTour);

void removeEdges(std::vector<Edge*>* toRemove, bool * canContainTour);

void addEdge(Edge * e);

void addEdges(std::vector<Edge*>* toAdd);

/*!
*  \brief Add an edge in the graph.
*  \param source : index of source node.
*  \param dest : index of dest node.
*  \param weight : weight of the edge.
*/
virtual  Edge * addEdge(int source, int dest, double weight);


void print() const;


graphNodes* getNodes();

graphEdges* getEdges();

/*!
 * \brief Return the number of nodes in this graph.
 */
int getSize() const;

/*!
 *  \brief Print the graph in graphviz dot format into file nomFichier.
 * A view of the graph can the be created using graphviz (www.graphviz.org).
 * type "dot -Tformat -o outfile grapheInDotFormatFile
 */
virtual void printDot(string nomFichier) const;




Node* removeNode(int nodeIndex);

void forceEdges(std::vector<Edge*>* toForce, bool *canContainTour);
void unforceEdges(std::vector<Edge*>* toForce);

int getNbArcForce();

void writeAssignmentData(string fileName);

};
#endif
