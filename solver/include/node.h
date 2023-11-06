#ifndef NODE_H 
#define NODE_H
#include <set>
#include <list>
#include <definitions.h>





/*!
 * \file node.h
 * \brief Header for Node an Edge.
 * \author Pascal Benchimol
 */




class Edge;
class Node;

class BiEdge;


//typedef std::set<int> nodeVoisins;
typedef std::list<Edge*> nodeEdges;
typedef std::list<Edge*> graphEdges;

/*!
 * \class Node
 * \brief A node of a graph.
 */

class Node {
private:
int index;
//nodeVoisins voisins; /*!< Set of neighboors*/
nodeEdges edges;	/*!< Set of edges connected to this node*/

int nbForcedVoisins;
public:
Node();
Node(int index);
/*!
*  \brief Add an edge concerning this node.
*  \param e : The edge to add.
*
*  If edge e touches this node, e is added to edges and the new neighboor of this node is added to voisins.
*/
nodeEdges::iterator addEdge(Edge * e);

void removeEdge(nodeEdges::iterator position) ;

void addForcedVoisin();
void removeForcedVoisin();

int getNbForcedVoisins();

int getIndex() const;

/*!
*  \brief Return set of neighboors of this node.
*
*/
//const nodeVoisins * getVoisins();


/*!
*  \brief Return set of edges connected to this node.
*
*/
nodeEdges * getEdges();





 /*!
 *  \brief Comparison operator.
 *
 *  a<b if
 *  a.index < b.index
 */
bool  operator<(const Node & b) const; 
/*!
*  \brief Comparison operator.
*
*  a==b if
*  a.index == b.index
*/
bool  operator==(const Node & b) const; 
bool  operator!=(const Node & b) const; 

};


/*!
 * \class Edge
 * \brief An edge of a graph.
 */
class Edge{
private:
Node* source; /*!< Source node*/
Node* dest; /*!< Destination node*/
double weight; /*!< Weight of this edge*/

bool removed;

double replacement_cost;

bool forced;

nodeEdges::iterator sourcePosition;
nodeEdges::iterator destPosition;
graphEdges::iterator graphPosition;

BiEdge * biedge;

public:

Edge(Node * source,  Node * dest, double weight);
Edge(Edge& e);

double coutMarg;

void force(bool * canContainTour);
void unforce();

bool isForced() const;

void setPositions(nodeEdges::iterator sourcePos, nodeEdges::iterator destPos,graphEdges::iterator graphPos);

nodeEdges::iterator getSourcePosition() const;
nodeEdges::iterator getDestPosition() const;
graphEdges::iterator getGraphPosition() const;

void remove();
void putBack();

bool isRemoved();


bool canBeRemoved();

double getReplacementCost();

void setReplacementCost(double cost);

/*!
*  \brief Return source node of this edge.
*
*/
Node * getSource() const;
/*!
*  \brief Return destination node of this edge.
*
*/
Node * getDest() const;
/*!
*  \brief Return weight of this edge.
*
*/
double getWeight();

double setWeight(double w);

 /*!
 *  \brief Comparison operator.
 *
 *  a<b if
 *  a.source < b.source
 * or a.source==b.source && a.dest < b.dest
 */
bool  operator<(const Edge & b) const; 
/*!
*  \brief Comparison operator.
*
*  a == b if
*  a.source==b.source && a.dest==b.dest
*/
bool  operator==(const Edge & b) const; 

BiEdge* getBiEdge();

void setBiEdge( BiEdge* _biedge);

void print() const;
};



#endif
