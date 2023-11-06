#include <stsp.h>
#include <limits>

SymTSP::SymTSP(int n, bool fromATSP) :
	Graphe::Graphe(n),
	isFromATSP(fromATSP)
{	
}


SymTSP::~SymTSP() 
{
}

void SymTSP::addEdge(Edge * e) {
	Graphe::addEdge(e);
}


Edge *  SymTSP::addEdge(int source, int dest, double weight) {

	if(source == dest) {
		cout << "Cannot add an edge from " << source << " to " << dest << endl;
		exit(0);
	}	

	if(source > dest) {
		Edge * e=Graphe::addEdge(dest,source,weight);
		
		return e;
	}
	else {
		Edge *  e=Graphe::addEdge(source,dest,weight);
		return e;
	}
	


}

Edge * SymTSP::addEdge(int source, int dest) {
	return addEdge(source, dest, 0.0);
}





void SymTSP::printDot(string nomFichier) const {
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"graph g{"<<endl;
    
      
    	
		for( graphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
			//Edge ed = *it;
			fb << (*it)->getSource()->getIndex() <<" -- " << (*it)->getDest()->getIndex()<<
				"[label = \" "<< (*it)->coutMarg<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		fb << "}"<<endl;


        // Ferme le fichier :
        fb.close();
    }


}


/*!
 * \class weightedNode
 * \brief Wrap a Node and store some data used in SymTSP::Prim() algorithm.
 *
 * Only used in SymTSP::Prim().
 */
class weightedNode {
public :
	Node* nodeRef;
    
	double weight;
	int parent;
	bool seen;
	
	weightedNode() :
	 weight(std::numeric_limits<double>::infinity()), parent(-1), seen(false)
	{};

	weightedNode(Node* ref, double w) :
		nodeRef(ref), weight(w), parent(-1), seen(false)
	{};
	    
	bool  operator<(const weightedNode & b) const {
		return this->weight > b.weight;    
	}
	
};

struct wnodeCompare {
 bool operator() (weightedNode* a, weightedNode* b) const
  {return *a<*b;}

};




/*!
* \class repNode
* \brief Object wrapping a Node. Only used in SymTSP::Kruskal() algorithm. Implements find-union mechanism for quick loop detection in Kruskal algorithm.
*
*/
class repNode{
	Node * nodeRef; /*!< Pointer to a Node*/
	int rank; /*!< Height of the cluster tree below this repNode*/
	repNode* rep; /*!< Pointer to the representative element of the clusted this repNode belongs to*/
	
	public :
	
	/*!
	 *  \brief Constructeur.
	 */
	repNode() :
		rank(0)
	 {rep=this;};
	
	/*!
	 * \brief Constructeur.
	 * \param n : pointeur to the node this repNode will wrap.
	 */
	repNode(Node* n) :
		nodeRef(n),
		rank(0)
	{rep=this;};
	
	
	/*!
	 *  \brief Find representative element of the cluster this repNode belongs to. Store the result in rep.
	 */
	repNode* findRep() {
		if(rep->nodeRef->getIndex()==nodeRef->getIndex()) {
			return rep;
		}
		rep = rep->findRep();
		return rep;
	};
	
	/*!
	 *  \brief Merge the cluster represented by this repNode with the cluster represented by b.
	 */
	void Union(repNode* b) {
		repNode* x = findRep();
		repNode* y = b->findRep();
		if(x->rank == y->rank) {
			x->rep = y->rep;
			y->rank=y->rank+1;
			return;
		}
		
		if(x->rank < y->rank) {
			x->rep=y->rep;
			return;
		}
		
		// else (x->rank > y->rank
		y->rep=x->rep;
		return;
		
	}

};	

//On classe les arcs par poids décroissant pour pouvoir accéder au plus petit
//élement par un appel à pop_back
bool compareEdgewithWeight (Edge* i,Edge* j) { 
	return (i->getWeight() > j->getWeight()); 
}

/*
SymTSP * SymTSP::Kruskal() {

	//On crée les repNodes et une liste pour y accéder.
	repNode* repNodes[maxSize]; 
	//Le tableau est de taille maxSize et pas getSize() pour pouvoir garder la possibilité d'enlever des noeuds au graphe
	
	for(std::map<int,Node>::iterator it = nodes.begin(); it!= nodes.end(); ++it) {	
		repNodes[it->first] = new repNode(&(it->second));
	}

	//On trie les arcs 
	std::vector<Edge*> sortedEdges;
	for(graphEdges::iterator it = edges.begin(); it!=edges.end(); ++it) {
		Edge * e =new Edge(**it);
		sortedEdges.push_back(e);
	}
	sort(sortedEdges.begin(), sortedEdges.end() , compareEdgewithWeight);
	
	
	SymTSP * arbre = new SymTSP(getSize(),false);
	double totalWeight=0.0;
	while(sortedEdges.size() > 0 ) {
		Edge * e = sortedEdges.back();
		//Les arcs sont triés par poids décroissant, 
		//donc pop_back donne le min
		sortedEdges.pop_back();
		
		repNode* a = repNodes[e->getSource()->getIndex()];
		repNode* b = repNodes[e->getDest()->getIndex()];
		
		//Si ca ne crée pas de cycle, on ajoute l'arc à l'arbre
		if(a->findRep() != b->findRep()) {
			a->Union(b);
			arbre->addEdge(e->getSource()->getIndex(), e->getDest()->getIndex(), e->getWeight());
			totalWeight+= e->getWeight();
		}
	
	}


	cout <<"Cost of the minimum spanning tree "<<totalWeight<<endl;	
	return arbre;

}
*/
/*
SymTSP * SymTSP::PrimMieux() {


		

	bool debug = false;
	
		
	Pointeur* added[getSize()]; //Ens des noeuds de l arbre
	
	
	
	std::vector<Pointeur*> toAdd; //tableau servant pour le heap
	
	
	
	toAdd.reserve(getSize()); //On connait deja la taille du vecteur : on reserve la mémoire
	
	Pointeur* pointeurs[maxSize]; // tableau pour acceder au wnodes
	//Le tableau est de taille maxSize et pas getSize() pour pouvoir garder la possibilité d'enlever des noeuds au graphe
	
	
	
	//initialisation

	//Fabrication des pointeurs
	int i=0;
	for(std::map<int,Node>::iterator it = nodes.begin(); it!= nodes.end(); ++it) {
		
		//Attention à bien initialiser les pointeurs avec leur position dans le heap
		Pointeur * p = new Pointeur(&(it->second), i, std::numeric_limits<double>::infinity());
		i++;
		
		pointeurs[it->first]=p;
		toAdd.push_back(p);
		
	}
	
	//fabrication du heap
	Heap * heap = new Heap(toAdd);
	heap->make_heap();
	
	if(debug) {
	  	cout << "heap contains: ";
		heap->print();
		cout << endl;
	}
	
	int lastAddIndex=0;
	while( heap->getTab()->size() >0) {
		
		if(debug) {
		  	cout << "heap contains: ";
			heap->print();
		}
		
		//On recupere l'element min du heap.
		Pointeur* p = heap->pop_heap();
		
		if(debug) {	cout<<"analyze "<<p->nodeRef->getIndex()<<endl;	}
		
		
		//on marque l'élement comme vu
		added[lastAddIndex]=p;
		lastAddIndex++;
		p->seen=true;
		int count=0;
		
		//Parcourir les voisins
		 for(     nodeEdges::const_iterator it(p->nodeRef->getEdges()->begin()),
	     itstop(p->nodeRef->getEdges()->end()); it!=itstop; ++it) {
	    	count++;
	    	int sourceNodeIndex = (*it)->getSource()->getIndex();
	    	int destNodeIndex = (*it)->getDest()->getIndex();
	    	int vref;
	    	//it->print();
	    	if(sourceNodeIndex==p->nodeRef->getIndex()) {
	    		vref=destNodeIndex;
	    	}
	    	else if (destNodeIndex==p->nodeRef->getIndex()) {
	    		vref=sourceNodeIndex;
	    	}
	    	else {
	    		cout << "Error in edges of node "<<p->nodeRef->getIndex() <<endl
	    			<<"Found an edge from "<<sourceNodeIndex <<" to "<<destNodeIndex<<endl;
    			exit(1);
	    	}
	    	
	    	//Attention : il faudrai vérifier que ce voisin est bien un 
	    	//element du graphe, mais prendrai log(nodes.size()) opérations...
	    	Pointeur * voisinPointeur = pointeurs[vref];
	    	
	    	double weight = (*it)->getWeight();
	    	//Si le voisin est plus proche
	    	if( !voisinPointeur->seen  &&  voisinPointeur->weight > weight) {
	    		//on change sa distance
		    	voisinPointeur->parent = p;
		    	voisinPointeur->weight=weight;
		    	//on met à jour le heap
		    	
		    //	cout <<"maj "<<voisinPointeur->nodeRef->getIndex()<< " old position "<<voisinPointeur->tabIndex;
		    	heap->monterDsHeap(voisinPointeur);
		    //	cout<<" new position "<<voisinPointeur->tabIndex<<endl;
			}
	    	
    	}
    	//cout <<"nb de voisins vu "<<count<<endl;
	
	}
	
	
	//fabrication de l'arbre résultant.
	SymTSP * arbre = new SymTSP(getSize(),false);

	double totalWeight =0.0;
	for(int i=1; i< getSize(); i++) {
		
		int source = (*(added+i))->nodeRef->getIndex();
		int dest = (*(added+i))->parent->nodeRef->getIndex();
		double weight =(*(added+i))->weight;
		totalWeight+= weight;
		if(debug) {
			cout << "node "<<source<<" has parent "<<dest<<endl;
		}
		arbre->addEdge(source, dest, weight);
	}
	

	cout <<"Cost of the minimum spanning tree "<<totalWeight<<endl;	
	return arbre;

}
*/
void SymTSP::printInTSPLIBFormat(string nomFichier) const {

	int n= getMaxSize();
	double matrix [n][n];
	
	for(int i=0;i <n; i++) {
		for(int j=0; j<n; j++) {
			matrix[i][j]=NO_EDGE;
		}
	}
	for(graphEdges::const_iterator it = edges.begin(); it!= edges.end(); ++it) {
		int s = (*it)->getSource()->getIndex();
		int t = (*it)->getDest()->getIndex();
		double w = (*it)->getWeight();
		
		if( s >t) {
			int ss=s;
			s=t;
			t=ss;
		}
		matrix[s][t]=w;
	}

	ofstream fb(nomFichier.c_str());
	// Teste si le fichier est ouvert :
	if (fb.is_open())
	{
		fb<<"TYPE: TSP"<<endl
			<<"DIMENSION: "<<n<<endl
			<<"EDGE_WEIGHT_TYPE: EXPLICIT"<<endl
			<<"EDGE_WEIGHT_FORMAT: UPPER_ROW"<<endl
			<<"EDGE_WEIGHT_SECTION"<<endl;
			
		for(int i=0; i< n; i++) {
			for(int j=i+1; j < n; j++) {
				fb<< matrix[i][j] << "  ";	
			}
			fb<<endl;
		}
		fb.close();
	
	}
}



