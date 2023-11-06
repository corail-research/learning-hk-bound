#ifndef BIGRAPHE_H
#define BIGRAPHE_H


#include<stsp.h>

#include <list>

#include <vector>
#include <binode.h>
#include <node.h>

#include <iostream>
#include <string>
#include <fstream>
//#include <stsp.h>
/*!
 * \file Graphe.h
 * \brief Header for Graphe.
 * \author Pascal Benchimol
 */

typedef std::vector<BiNode*> bigraphNodes;


using namespace std;

/*!
 * \class Graphe
 */
class BiGraphe {
private:
int size;
int matchingSize;
bigraphNodes rightNodes; /*!< Vector of nodes*/
bigraphNodes leftNodes; /*!< Vector of nodes*/
bigraphEdges edges;	/*!< List of edges*/


BiEdge* addEdge(int leftIndex, int rightIndex, double weight, Edge* _edge);



public:

/*!
 *  \brief Constructeur
 *  \param size : number of nodes.
 */
BiGraphe(SymTSP * stsp);

~BiGraphe();


/*!
 * \brief Return the maximum index of a node that can be found in this graphe.
 */
int getSize() const;
int getMatchingSize() const;


bigraphNodes* getRightNodes();
bigraphNodes* getLeftNodes();
bigraphEdges* getEdges();




void print() const;


void addEdge(BiEdge * edge);

void removeEdge(BiEdge * e);

/*!
 *  \brief Print the graph in graphviz dot format into file nomFichier.
 * A view of the graph can the be created using graphviz (www.graphviz.org).
 * type "dot -Tformat -o outfile grapheInDotFormatFile
 */
void printDot(string nomFichier) const;








};
#endif
