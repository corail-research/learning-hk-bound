#ifndef BINODE_H 
#define BINODE_H
#include <graphe.h>
#include <set>
#include <list>
#include <limits>
#include <definitions.h>
#include <string>
#include <iostream>
#include <sstream>
#include <node.h>




/*!
 * \file node.h
 * \brief Header for Node an Edge.
 * \author Pascal Benchimol
 */




class BiEdge;
class BiNode;


//typedef std::set<int> nodeVoisins;
typedef std::list<BiEdge*> binodeEdges;
typedef std::list<BiEdge*> bigraphEdges;

/*!
 * \class Node
 * \brief A node of a graph.
 */

class BiNode {
private:
const int index;
Node  * const   masterNode;
const bool set;   //true if in left set X, false if in right set Y
binodeEdges * edges;	/*!< Set of edges connected to this node*/



double potential;
BiEdge * parentEdge; //edge leading to parent in alternate tree
BiEdge* matchingEdge;

double slack;
binodeEdges minSlackEdges;

bool inTree;

list<BiNode*>::iterator position_listNotInTree;




public:


BiNode(int _index, bool _isInLeftSet,double _potential, Node * masterNode);
~BiNode();

int getIndex() const;
bool isInLeftSet() const;

Node * getMasterNode();

binodeEdges* getEdges();
binodeEdges::iterator addEdge(BiEdge * e);
void removeEdge(binodeEdges::iterator position) ;

double getPotential() const;
void setPotential(double _potential);

double getSlack() const;
void setSlack(double _slack);

void addMinSlackEdge(BiEdge* edge);
binodeEdges* getMinSlackEdges();

void init();

BiEdge * getParentEdge();
void setParentEdge(BiEdge * node);


BiEdge* getMatchingEdge();
void setMatchingEdge(BiEdge * edge);

void clear();

bool isInTree();
void putInTree();
void removeFromTree();



list<BiNode*>::iterator getPosition();
void setPosition(list<BiNode*>::iterator pos);

/*!
*  \brief Return set of edges connected to this node.
*
*/


string print() const;

};


/*!
 * \class Edge
 * \brief An edge of a graph.
 */
class BiEdge{
private:
BiNode* leftNode; /*!< Source node*/
BiNode* rightNode; /*!< Destination node*/
const double weight; /*!< Weight of this edge*/

Edge* edge;

bool removed;

binodeEdges::iterator leftPosition;
binodeEdges::iterator rightPosition;
bigraphEdges::iterator graphPosition;



public:
//Edge();
BiEdge(BiNode * _leftNode,  BiNode * _rightNode, double _weight, Edge* _edge);



bool isTight() const;





void setPositions(binodeEdges::iterator leftPos, binodeEdges::iterator rightPos, bigraphEdges::iterator graphPos);

binodeEdges::iterator getLeftPosition() const;
binodeEdges::iterator getRightPosition() const;
bigraphEdges::iterator getGraphPosition() const;


Edge* getEdge() const;



bool isRemoved();

void remove();

void putBack();

/*!
*  \brief Return source node of this edge.
*
*/
BiNode * getLeftNode() const;
BiNode * getRightNode() const;

/*!
*  \brief Return weight of this edge.
*
*/
double getWeight() const;

double getReducedCost() const;

std::string print() const;
};



#endif
