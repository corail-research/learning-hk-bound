#ifndef HEAP_H
#define HEAP_H
#include<vector>
#include <limits>
#include <math.h>
#include <node.h>





/*!
 * \class Pointeur
 * \brief Used in Heap. Wrap a Node and store its position in the heap.
 * Represent a node in the minimum spanning tree AND
 * the edge from this node to its parent in the minimum spanning tree.
 */
class Pointeur {
	public: 
	Node* nodeRef; /*!< Pointer to the Node this object wraps.*/

	int tabIndex; /*!< Index of this Pointeur in the heap's array. */

	double weight; /*<! Distance weight in Prim algorithm. */
	double realWeight;
	Pointeur* parent; /*<! Parent in the tree resulting from Prim algorithm. */
	bool seen; /*! True if this node is already in the tree in Prim algorithm. */
	
	int rank; /*<! Depth of this node in the minimum spanning tree */
	int degre; /*<! Degree of this node in the minimum spanning tree */
	

	
	std::vector<Pointeur*> fils; /*<! List of sons in the tree */
	
	std::vector<Edge* > supportCycleEdges;
	Edge* minEdgeInSupport;
	double minEdgeWeight;



	Edge* edge;
	
	void clear() {
		tabIndex=-1;
		parent=NULL;
		weight = std::numeric_limits<double>::infinity();
		realWeight = std::numeric_limits<double>::infinity();
		rank = -1;
		degre=0;
		seen=false;
		fils.clear();
		edge=NULL;
		supportCycleEdges.clear();
		minEdgeInSupport=NULL;
		minEdgeWeight=std::numeric_limits<double>::infinity();
	};

	/**
	* \brief Standard constructor.
	*/
	Pointeur() :
	tabIndex(-1),
	weight(std::numeric_limits<double>::infinity()),
	realWeight(std::numeric_limits<double>::infinity()),
	seen(false),
	parent(NULL),
	edge(NULL),
	minEdgeInSupport(NULL),
	minEdgeWeight(std::numeric_limits<double>::infinity())
	{};
	
	/**
	* \brief Standard constructor.
	* \param nr Pointer to the Node this object wraps.
	* \param tIndex Index of this Pointeur in the heap's array.
	* \param w Distance weight in Prim algorithm. Usually initialized at infinity.
	*/
	Pointeur(Node * nr, int tIndex, double w) :
	nodeRef(nr),
	tabIndex(tIndex),
	weight(w),
	realWeight(w),
	parent(NULL),
	seen(false),
	edge(NULL)
	{};
	
	bool  operator<(const Pointeur & b) const {
	
		if(weight < b.weight)
			return true;
		
/*
		  srand ( time(NULL) );
		int iSecret = rand() % 100;
		
		if(weight==b.weight &&  iSecret <50) {
			return true;
		}
*/
			
		/*
		if(weight==b.weight && degre < b.degre) {
			return true;
		}
		*/
		/*
		if(weight==b.weight && degre == b.degre && nodeRef->getIndex() <b.nodeRef->getIndex() ) {
			return true;
		}
		*/
		return false;
	
	}
	
	void addEdgeInSupport(Edge* ed, double we) {
		supportCycleEdges.push_back(ed );
		if(we < minEdgeWeight) {
			minEdgeWeight=we;
			minEdgeInSupport=ed;
		}
	}

};

/*!
 * \class Heap
 * \brief A custom made heap which allows decrease_key operation in O(log n).
 *
 *
 * This heap store Nodes wrapped in a Pointeur object. <br>
 * Every time a Node is moved in the heap, its Pointeur record its position in the heap.
 * Knowing the position of a Node allows decrease_key operation by moving up this Node in the heap
 * from its position. <br>
 *<br>
 * This heap use an array of Pointeur* to store its data.
 */
class Heap {

	private :
	std::vector< Pointeur* > tableau; /*<! Array containing data for the heap. */
	
	int getParent(int i);
	
	
	
	/**
	*	\brief Move down the pointeur p until it reaches its right position in the heap.
	* \param p Pointeur to move down.
	*/
	void descendreDsHeap(Pointeur* p);
	
	
	
	public :
	/**
	* \brief Standard constructor.
	*
	* Warning : pointeur->tabIndex must be initialized as required before calling this constructor.
	* \param tab Vector of pointeur to use in the heap.
	*/
	Heap(std::vector<Pointeur*> & tab);
	~Heap();
	
	
	/**
	*\brief Return the array containing data of this heap.
	*/
	std::vector< Pointeur*> * getTab() ;
	
	
	/**
	*\brief Order data in the array so as to it respects heap properties.
	*
	*Simply call descendreDsHeap() from n/2 to 0. <br>
	*<br>
	*Expected complexity : O(n).
	*/
	void make_heap();
	
	/**
	*\brief Return minimum element of the heap and order the remaining data so as it respects heap properties.
	*
	*Remove the root of the heap. Replace it with the last element in the heap, the call descendreDsHeap() on it.
	<br>
	*<br>
	*Expected complexity : O(log(n)).
	*/
	Pointeur * pop_heap();
	
	/**
	*\brief Add pointeur p to the heap and order the data so as it respects heap properties.
	* p is appended at the end of the heap (well the array) and monterDsHeap() is called on it.
	*<br><br>
	*Expected complexity : O(log(n))
	* \param p Pointeur to add into the heap.
	*/
	void push_heap( Pointeur * p);
	
	/**
	*\brief Move up pointeur p until it reaches its right place in the heap.

	* \param p Pointeur to move up.
	*/
	void monterDsHeap(Pointeur * p);
	
	void print();

	
};




#endif
