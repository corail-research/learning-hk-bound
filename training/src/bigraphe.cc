#include <bigraphe.h> 
#include <vector>


BiGraphe::BiGraphe(SymTSP* stsp) 

{


	int n=stsp->getMaxSize();
	if(!stsp->isFromATSP) {
		 return;
		size=n;
	 	matchingSize=n;
	
	
		//Node creation
		leftNodes.reserve(n);
		rightNodes.reserve(n);
		for(graphNodes::iterator it = stsp->getNodes()->begin(); it != stsp->getNodes()->end(); ++it) {
			int index = (*it)->getIndex();
			
			
			
			leftNodes.push_back( new BiNode(index,true,0, (*it)));
	
			rightNodes.push_back( new BiNode(index,false,0, (*it)));
			
		}
	
	
		//cout<<"Nodes created"<<endl;
		//Edges creation
		for(graphEdges::iterator it = stsp->getEdges()->begin(); it != stsp->getEdges()->end(); ++it ) {
			const int sourceIndex= (*it)->getSource()->getIndex();
			const int destIndex = (*it)->getDest()->getIndex();
			const double weight= (*it)->coutMarg;
			//const double weight= (*it)->getWeight();
		
			//We add two biEdge for each edge
			
			addEdge(sourceIndex, destIndex, weight, (*it));
		
			addEdge(destIndex, sourceIndex, weight, (*it));
		
		}
	
		//cout<<"Edges created"<<endl;
		
	
	}
	else {
		cout<<"ATSP graph"<<endl;
		if(n %2 != 0) {
			cout<<"Error, graphe is from ATSP and number of nodes is not even"<<endl;
			exit(1);
		}
		
		n=n/2;
		size=n;
		matchingSize=n;
		//Node creation
		leftNodes.reserve(n);
		rightNodes.reserve(n);
		
		for(graphNodes::iterator it = stsp->getNodes()->begin(); it != stsp->getNodes()->end(); ++it) {
			int index = (*it)->getIndex();
			
			//edges in the stsp built from ATSP
			//in nodes : i
			//out nodes : size +i
			
			//left nodes are for outgoing edges
			if(index >= n) {
				leftNodes.push_back( new BiNode(index-n,true,0, (*it)));
			}
			//right nodes are for incoming edges
			else {
				rightNodes.push_back( new BiNode(index,false,0, (*it)));
			}
		}
	
		
		//cout<<"Nodes created"<<endl;
		//Edges creation v1
		
		for(graphEdges::iterator it = stsp->getEdges()->begin(); it != stsp->getEdges()->end(); ++it ) {
			int sourceIndex= (*it)->getSource()->getIndex();
			int destIndex = (*it)->getDest()->getIndex();
			
			if(sourceIndex > destIndex) {
				int ss = sourceIndex;
				sourceIndex=destIndex;
				destIndex=ss;
			}
			
			//we skip edges added during atsp -> stsp transformation
			if(destIndex==sourceIndex+n) {
				continue;
			}
			destIndex-=n;
			//double weight= (*it)->coutMarg;
			double weight= (*it)->getWeight();
		
			//We add one biEdge for each edge
			
			//edges in the stsp built from ATSP
			//in nodes : i
			//out nodes : size +i
			//so edges in ATSP are from destIndex to sourceIndex since destIndex > sourceIndex
			
			BiEdge * biedge =addEdge(destIndex, sourceIndex, weight, (*it));
			(*it)->setBiEdge(biedge);
			
		
		
		}
		
		//cout<<"Edges created"<<endl;
		
	}

	
}




BiGraphe::~BiGraphe() {
	for(bigraphNodes::iterator it = leftNodes.begin(); it!= leftNodes.end(); ++it) {
		delete (*it);
	
	}
//	delete &Snodes;


	for(bigraphNodes::iterator it = rightNodes.begin(); it!= rightNodes.end(); ++it) {
		delete (*it);
	
	}
//	delete &Tnodes;
	
	
	for(bigraphEdges::iterator it = edges.begin(); it!= edges.end(); ++it) {
		delete (*it);
	
	}
//	delete &edges;
}

int BiGraphe::getSize() const {
	return size;
}


int BiGraphe::getMatchingSize() const {
	return matchingSize;
}





BiEdge * BiGraphe::addEdge(int leftIndex, int rightIndex, double weight, Edge* _edge) {

	/*
	BiNode * leftNode = NULL;
	for(bigraphNodes::iterator it = leftNodes.begin(); it!= leftNodes.end(); ++it) {
		if( (*it)->getIndex() == leftIndex) {
			leftNode = (*it);
			break;
		}
	}
	if(leftNode == NULL) {
		cout<<"Error while creating bipartite graphe. Cannot find left node "<<leftIndex<<endl;
		exit(1);
	}
	
	BiNode * rightNode = NULL;
	
	for(bigraphNodes::iterator it = rightNodes.begin(); it!= rightNodes.end(); ++it) {
		if( (*it)->getIndex() == rightIndex) {
			rightNode = (*it);
			break;
		}
	}
	if(leftNode == NULL) {
		cout<<"Error while creating bipartite graphe. Cannot find right node "<<rightIndex<<endl;
		exit(1);
	}
	*/

	BiNode * leftNode = leftNodes.at(leftIndex);
	BiNode * rightNode = rightNodes.at(rightIndex);

	BiEdge * edge = new BiEdge(leftNode, rightNode, weight, _edge);
	
	
	addEdge(edge);

	return edge;

}

void BiGraphe::addEdge(BiEdge * edge) {
	if(edge->isRemoved() == false) {
		return;
	}
	
	
	//cout<<"Add edge "<<edge->print()<<endl;
	edges.push_back(edge);
	bigraphEdges::iterator graphePos = --edges.end();
	binodeEdges::iterator leftPos   = edge->getLeftNode()->addEdge(edge);
	binodeEdges::iterator rightPos   = edge->getRightNode()->addEdge(edge);

	edge->setPositions(leftPos,rightPos, graphePos);
	
	edge->putBack();
}


void BiGraphe::removeEdge(BiEdge * e) {
	if(e->isRemoved() == true) {
		return;
	}
	//cout<<"Remove edge "<<e->print()<<" "<<e->getEdge()->isForced()<<endl;
	e->getLeftNode()->removeEdge(e->getLeftPosition());
	e->getRightNode()->removeEdge(e->getRightPosition());
	edges.erase(e->getGraphPosition());
	e->remove();
}



bigraphNodes* BiGraphe::getLeftNodes() {
	return &leftNodes;
}
bigraphNodes* BiGraphe::getRightNodes() {
	return &rightNodes;
}


bigraphEdges* BiGraphe::getEdges() {
	return &edges;
}






void BiGraphe::print() const{
	
	cout << "Nb of nodes : "<< size << endl;
	cout << "Nb of edges : "<< edges.size() << endl;
	for( bigraphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
		cout<<(*it)->print()<<endl;
	}
}





void BiGraphe::printDot(string nomFichier) const {
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"digraph g{"<<endl;
    	fb<<"rankdir=LR"<<endl;;
      
    	
		for( bigraphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
			//Edge ed = *it;

		//	if( (*it)->getSourceNode()->isInLeftSet()) {
		//		continue;
		//	}
			
			fb << (*it)->getLeftNode()->print() <<" -> " << (*it)->getRightNode()->print()<<
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






