#include <graphe.h> 
#include <vector>

bool printMST=false;
static bool debug=false;

extern bool verbose;

Graphe::Graphe(int s)  :
 maxSize(s)
{
	
	nodes.reserve(s);
	for(int i=0; i<s; i++) {

		nodes.push_back (new Node(i) );
	}
}


Graphe::Graphe()
	{};

Graphe::~Graphe() {
	
	for(graphNodes::iterator it = nodes.begin(); it != nodes.end(); ++it) {
		delete (*it);
	}
	

	for(graphEdges::iterator it = edges.begin(); it!= edges.end(); ++it) {
		delete *it;
	
	}
	
}

int Graphe::getMaxSize() const {
	int n=maxSize;
	return n;
}




void Graphe::removeEdge(Edge * e, bool * canContainTour) {
	//check si on peut enlever cet arc
	
	if( e->isRemoved())  {
		if(verbose) cout<<"Error, try to remove a removed edge...";
		if(verbose) e->print();
		if(verbose) cout<<endl;
		return;
	}
	
	
	if(e->getSource()->getEdges()->size() <= 2) {
		//cout<<"Cannot remove edge from "<<e->getSource()->getIndex()<<" to "<<e->getDest()->getIndex()
		//<<" because node "<<e->getSource()->getIndex()<<" has only "<<e->getSource()->getEdges()->size()<<" edges."<<endl;
		//exit(1);
		*canContainTour = false;
	}
	if(e->getDest()->getEdges()->size() <= 2) {
		//cout<<"Cannot remove edge from "<<e->getSource()->getIndex()<<" to "<<e->getDest()->getIndex()
		//<<" because node "<<e->getDest()->getIndex()<<" has only "<<e->getDest()->getEdges()->size()<<" edges."<<endl;
		//exit(1);
		*canContainTour = false;
	}

	if(debug) {
		if(verbose) cout<<"Remove edge from "<<e->getSource()->getIndex()<<" to "<<e->getDest()->getIndex()<<endl;
		
	}
	e->getSource()->removeEdge(e->getSourcePosition());
	
	e->getDest()->removeEdge(e->getDestPosition());
	edges.erase(e->getGraphPosition());
	
	
	e->remove();
	
	//delete e;
}

void Graphe::removeEdges(std::vector<Edge*>* toRemove, bool * canContainTour) {
	* canContainTour = true;
	for(std::vector<Edge*>::iterator it = toRemove->begin(); it != toRemove->end(); it++) {
		removeEdge(*it, canContainTour);
		
	}
	if(verbose) cout<<toRemove->size()<<" edges removed. Remaining "<<edges.size()<<endl;
}

Edge * Graphe::addEdge(int source, int dest, double weight) {

    	// Check if an edge can be added from source to dest.
	if( source >= getSize() || source <0) {
		cout << "asked for source node number "<< source << " in a graph of dimension "<<getSize()<<endl;
		exit(1);
	}
	if( dest >= getSize() || dest <0) {
		cout << "asked for dest node number "<< dest << " in a graph of dimension "<<getSize()<<endl;
		exit(1);
	}

	if(source == dest) {
		cout << "Cannot add an edge from " << source << " to " << dest << endl;
		exit(1);
	}	

	Node* sourceNode = nodes.at(source);
	Node* destNode = nodes.at(dest);
	
	Edge * e = new Edge(sourceNode, destNode, weight);


	/*
	//Check if an edge from source to dest already exists.

	std::pair<graphEdges::const_iterator, bool> exist= edges.insert(e);
	if(!exist.second) {
		if(exist.first->getWeight()!=weight) {
			cout << "Warning : Edge from " << source << " to "<< dest<< " already exist"<<endl
			<<"stored value : " << exist.first->getWeight() <<". Tried to add value "<<weight<< endl;
		}
	}
	*/
	addEdge(e);
	
	return e;

}

void Graphe::addEdge(Edge* e) {
	
	/*
	if( containsEdge(e->getSource()->getIndex()
	, e->getDest()->getIndex() ) ) {
		cout<<"Edge added 2 times :"<<endl;
		e->print();
		exit(1);
	}
	*/
	//Insert new edge in edges
	if(e->isRemoved() == false ) {
		return;
		//cout<<"Error. Adding an edge which is already in the graph."<<endl;
		//e->print();
		//exit(1);
	}
	
	edges.push_back(e);
	e->putBack();

	graphEdges::iterator graphPos = --edges.end();
	
	
	//Add edge to the nodes.
	nodeEdges::iterator sourcePos = e->getSource()->addEdge(e);
	nodeEdges::iterator destPos = e->getDest()->addEdge(e);
	e->setPositions(sourcePos,destPos, graphPos);

}

void Graphe::addEdges(std::vector<Edge*> *toAdd) {
	for(std::vector<Edge*>::iterator it = toAdd->begin(); it != toAdd->end(); it++) {
		addEdge(*it);
	}

}


graphNodes* Graphe::getNodes() {
	return &nodes;
}

graphEdges* Graphe::getEdges() {
	return &edges;
}

int Graphe::getSize() const{
	return nodes.size();
}


/*
Node* Graphe::removeNode(int indexNode) {
	map<int,Node>::iterator result = nodes.find(indexNode);
	if(result==nodes.end()) {
		cout <<"Node "<<indexNode<<" is not in the graph. Cannot remove it."<<endl;
		exit(1);
	}
	
	nodes.erase(result);
	
	
	for(set<Edge>::iterator it = result->second.getEdges()->begin();
						it !=result->second.getEdges()->end(); ++it) {
		
		int erased=edges.erase(*it);					
		if(erased==0) {
			cout <<"Node "<<indexNode<<" is not in the graph. Cannot remove it."<<endl;
			exit(1);
		}
	}
}
*/
void Graphe::print() const{
	
	cout << "Nb of nodes : "<< getSize() << endl;
	for( graphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
		Edge ed = **it;
		ed.print();
	}
}




void Graphe::printDot(string nomFichier) const {
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"graph g{"<<endl;
    
      
    	
		for( graphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
			//Edge ed = *it;
			fb << (*it)->getSource()->getIndex() <<" -- " << (*it)->getDest()->getIndex()<<
				"[label = \" "<< (*it)->getWeight()<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		fb << "}"<<endl;


        // Ferme le fichier :
        fb.close();
    }


	


}

void Graphe::forceEdges(std::vector<Edge*>* toForce, bool *canContainTour) {
	int nbForce=0;

	for(std::vector<Edge*>::iterator it = toForce->begin(); it!= toForce->end(); ++it) {
		if( ! (*it)->isForced() ) {
			nbForce++;
			(*it)->force(canContainTour);
			if(debug)
				cout<<"Force edge from "<<(*it)->getSource()->getIndex()<<" to "<<(*it)->getDest()->getIndex()<<endl;
		}
	}
	if(verbose) cout<<nbForce<<" edges forced."<<endl;

}


void Graphe::unforceEdges(std::vector<Edge*>* toForce) {
	for(std::vector<Edge*>::iterator it = toForce->begin(); it!= toForce->end(); ++it) {
		(*it)->unforce();
	}
}


int Graphe::getNbArcForce() {
	int result=0;
	for(list<Edge*>::iterator it = edges.begin(); it!= edges.end(); ++it) {
		if((*it)->isForced())
			result++;
	}
	
	return result;
}


void Graphe::writeAssignmentData(string fileName) {
	  ofstream fb(fileName.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb<<"data;"<<endl;
    	
    	//write set of nodes
    	fb <<"set Nodes :=";
    	for(graphNodes::iterator it = nodes.begin(); it!= nodes.end(); ++it) {
    		fb<<(*it)->getIndex()<<" ";
    	
    	}
    	fb<<";"<<endl<<endl
    	//write set of edges
    	<<"set Edges :=";
    	for(graphEdges::iterator it = edges.begin(); it!= edges.end(); ++it) {
    		fb<<(*it)->getSource()->getIndex()<<"_"<<(*it)->getDest()->getIndex()<<" ";
    	}
    	fb<<";"<<endl<<endl;
    	
    	//write set of neighboors
    	for(graphNodes::iterator it = nodes.begin(); it!= nodes.end(); ++it) {
    		Node* n = (*it);
    		fb<<"set V["<<n->getIndex()<<"] :=";
    		
    		for(nodeEdges::const_iterator it = n->getEdges()->begin(); it != n->getEdges()->end();++it) {
    			fb<<(*it)->getSource()->getIndex()<<"_"<<(*it)->getDest()->getIndex()<<" ";
    		
    		}
    		fb<<";"<<endl;
    	
    	}
    	
    	fb<<endl;
    	
    	//write cost
    	fb<<"param cost :="<<endl;
    	for(graphEdges::iterator it = edges.begin(); it!= edges.end(); ++it) {
    		fb<<(*it)->getSource()->getIndex()<<"_"<<(*it)->getDest()->getIndex()<<" "<<(*it)->coutMarg<<endl;
    		//fb<<(*it)->getSource()->getIndex()<<"_"<<(*it)->getDest()->getIndex()<<" "<<(*it)->getWeight()<<endl;
    	}
    	fb<<";"<<endl
    	<<"end;"<<endl;
    	
    	
	}
}

