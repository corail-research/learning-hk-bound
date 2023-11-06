#ifndef STSP_H 
#define STSP_H
#include <graphe.h>
#include <heap.h>
#include <vector>
#include <algorithm>

//class Arbre;



	/*! \class SymTSP
	* \brief A Symmetric Travelling Salesman Problem object.
	*
	* Symmetry of these instances is used to reduce the number of stored edges.
	* The only edges stored are from i to j where i<j.
	* All methods with regards to edges are accordingly overrided.
	*/

class SymTSP : public Graphe {

public :
 /*!
 *  \brief Constructeur
 *  \param n : number of nodes.
 */
 
 const bool isFromATSP;
 
SymTSP(int n, bool fromATSP);
~SymTSP();
void addEdge(Edge * e);
Edge * addEdge(int source, int dest, double weight);
Edge * addEdge(int source, int dest);
bool containsEdge(int source, int dest) const;
const Edge* getEdge(int source, int dest) const;

void printDot(string nomFichier) const;

/*!
 * \brief Implements Prim algorithm using stl. <br>
 *
 *
 *Informations needed for each node (distance, parent) are wrapped into weightedNode objects. <br><br>
 * This pseudo-code is used: <br>
 * \verbatim
 * d[n] : distance of node n 
 * d[n] <- infinity 
 * p[n] : parent of node n in the minimum spanning tree 
 * p[n] <- null 
 * heap H contains all nodes, sorted according to distance d 
 * 
 * while( H is not empty) { 
 *		n = H.pop_min(); 
 *		forall(node v adjacent to n) { 
 *			if( v is in H && d[v] > d[n] + cost(n,v) ) { 
 *				d[v] <- d[n] + cost(n,v); 
 *				p[v] <- n; 
 *			} 
 *		} 
 *		if(d[v] has been updated for some v) { 
 *			make_heap(H);
 *		} 
 * } 
 * \endverbatim
 * Because there is no decrease_key_in_heap algorithm in stl and there is no way to find the position
 * of an element in the heap (without going throught the whole heap), the function make_heap has
 * to be used after modification of d[v]. <br> <br>
 * There is |V| iteration of while loop. <br>
 * In each iteration : <br>
 * pop_min takes O(log (|V|) ) <br>
 * there is on average |E|/|V| neighboors <br>
 * make_heap takes O(|V|) <br>
 * So it takes O( |V| ( log (|V|) + |E|/|V| + |V| ) = O( |E|+|V|² + |V|*log(|V|) ) operations <br> <br>
 *
 * The expected complexity is O( |E|+ |V|*log(|V|) )
 */
SymTSP * Prim();

/*!
 * \brief Implements Prim algorithm using a custom made Heap. <br>
 *
 *
 * The custom made Heap allows to decrease a key in O(log n), and so improve
 * complexity of Prim algorithm. <br>
 *
 
 *
 * *Informations needed for each node (distance, parent, position in heap) are wrapped into Pointeur objects. <br><br>
  
 * This pseudo-code is used: <br>
 * \verbatim
 * d[n] : distance of node n 
 * d[n] <- infinity 
 * p[n] : parent of node n in the minimum spanning tree 
 * p[n] <- null 
 * heap H contains all nodes, sorted according to distance d 
 * while( H is not empty) { 
 *		n = H.pop_min(); 
 *		forall(node v adjacent to n) { 
 *			if( v is in H && d[v] > d[n] + cost(n,v) ) { 
 *				d[v] <- d[n] + cost(n,v); 
 *				p[v] <- n; 
 *				H.decrease_key(v) 
 *			} 
 *		} 
 * } 
 * \endverbatim
 * There is |V| iteration of while loop. <br>
 * In each iteration : <br>
 * pop_min takes O(log (|V|) ) <br>
 * there is on average |E|/|V| neighboors <br>
 * each decrease_key takes takes O(log(|V|)) <br>
 * So it takes O( |V| ( log (|V|) + |E|/|V|*log(|V|)  ) ) = O( (|E|+ |V|) log (|V|) ) operations <br> <br>
 * The expected complexity is O( (|E|+|V|) * log(|V|) )
 */
SymTSP * PrimMieux();


/*!
 * \brief Implements Kruskal algorithm using a find-union mechanism to detect loops. <br>
 *
 *
  *Informations needed for each node (representative) are wrapped into repNode objects. <br><br>
 *
 * This pseudo-code is used: <br>
 * \verbatim
 * H : vector containing all edges sorted by edge weight 
 * rep[n] : a node representing the cluster n belongs to. 
 * rep[n] <- n 
 * while( H is not empty) { 
 *		a = e.source; 
 *		b = e.destination; 
 *		if(  FIND(a) != FIND(b) ) { 
 *			UNION(a,b); 
 *			add edge e to the tree 
 *		} 
 * } 
 * \endverbatim <br>
 * See class repNode for implementation of FIND and UNION.<br> <br>
 * Sorting of |E| edges is done in O(|E|log(|E|) according to stl documentation <br>
 * H contains |E| edges, so there will be |E| iteration of while loop <br>
 * In the while loop we do two calls to FIND  and one to UNION. <br>
 * We know that complexity of |E| calls to FIND or UNION has complexity  O( |E| alpha(|E|,|V|) ) <br> <br>
 * The expected complexity is O(|E|log(|E|+|E| alpha(|E|,|V|) )
 */
SymTSP * Kruskal();

void printInTSPLIBFormat(string nomFichier) const;

};
#endif
