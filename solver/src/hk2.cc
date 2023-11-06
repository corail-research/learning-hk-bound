#include <hk2.h>
extern  bool debug;






/*
* Return the list of edges of the MST which can't be removed.
* If one try to remove them, the HKBound will rise above the given upperBound.
*
*/
void HeldKarp::forceMSTEdges(double upperBound, std::vector<Edge*>* toForce) {
	


	//look at edges of the MST
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end(); ++it) {
		if((*it)->rank == 0 ){ //on ne regarde pas la racine de l'arbre
			continue;
		}
		
		if((*it)->nodeRef->getIndex() == index1node ){ //on ne regarde pas le oneNode
			continue;
		}
		
		if((*it)->edge== NULL ) {
			cout<<(*it)->nodeRef->getIndex()<<endl;
			cout<<(*it)->rank<<endl;
			cout<<"erreur, node sans parent dans le mst"<<endl;
			exit(1);
		}
		
		
		if((*it)->edge->isForced() ) {//on ne regarde pas les arcs deja forcé
			continue;
		}
		
		double coutMarg= (*it)->minEdgeWeight - (*it)->realWeight;
		
		(*it)->edge->setReplacementCost(coutMarg);
		
		if(coutMarg+HKBound - upperBound > PRECISION ) {
			toForce->push_back((*it)->edge);
			
		if(debug) {
			cout<<"Edge from "<<(*it)->nodeRef->getIndex()<<" to "<<(*it)->parent->nodeRef->getIndex()<<"  "<<(*it)->realWeight
				<<" can be forced : "
				<<"cout marg + HKBound = "<<(coutMarg+HKBound)<<" > UB"<<endl;
			}
		}	
	}
	
	//look at edges connected to the one-node
	if(oneEdge3Weight < std::numeric_limits<double>::infinity()) {
	
		if(!oneEdge1->isForced()) {
			double coutMarg=oneEdge3Weight - oneEdge1Weight;
			oneEdge1->setReplacementCost(coutMarg);
			if(coutMarg+HKBound - upperBound > PRECISION ) {
				toForce->push_back(oneEdge1);
			}

		}
		if(!oneEdge2->isForced()) {
			double coutMarg=oneEdge3Weight - oneEdge2Weight;
			oneEdge2->setReplacementCost(coutMarg);
			if(coutMarg+HKBound - upperBound > PRECISION ) {
				toForce->push_back(oneEdge2);
			}

		}
	}
	
//	cout<<result->size()<<" edges can be forced because of filtering"<<endl;
	return;

}






//If a node has degree 2, then its two edges can be forced.
void HeldKarp::forceEdgesWithDegree2(std::vector<Edge*>* toForce) {
	
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end() ; ++it) {
		if((*it)->nodeRef->getEdges()->size() != 2 ) {
			continue;
		}
		Node * node = (*it)->nodeRef;
		
		for(nodeEdges::const_iterator vit = node->getEdges()->begin(); vit != node->getEdges()->end(); ++vit ) {
			if( !(*vit)->isForced() ) {
				toForce->push_back((*vit));
				
				if(debug) {
					cout<<"Node "<<node->getIndex()<<" has degree 2, so edge ";
					(*vit)->print();
					cout<<" is forced. "<<endl;
				}
			
			}
			
		}

		
	}
	return;
}






void HeldKarp::PrimForce(bool * containsTour, std::vector<Edge*>* toForce) {

	bool debug = false;
	
	vector<Pointeur*> pointeurs2=vector<Pointeur*>(graphe->getMaxSize());
	
	for(graphNodes::iterator it = graphe->getNodes()->begin(); it!= graphe->getNodes()->end(); ++it) {
			pointeurs2.at((*it)->getIndex())=new Pointeur((*it), -1, std::numeric_limits<double>::infinity());
	}
	
	*containsTour = true;

	
	
	std::vector<Pointeur*> added; //Ens des noeuds de l arbre
	added.reserve(graphe->getSize()); 
	
	std::vector<Pointeur*> toAdd; //tableau servant pour le heap
	toAdd.reserve(graphe->getSize()); //On connait deja la taille max du vecteur : on reserve la mémoire
	
	//fabrication du heap
	Heap heap(toAdd);
	

	//initialize heap with a node 

	heap.push_heap( (*pointeurs2.begin()) );	

	int cutSize = 0;
	
	if(debug) {
	  	cout << "heap contains: ";
		heap.print();
		cout << endl;
	}



	while( true) {
		
		if(debug) {
		  	cout << "heap contains: ";
			heap.print();
		}
		
		//On recupere l'element min du heap.
		Pointeur* p = heap.pop_heap();

		//on marque l'élement comme vu
		added.push_back(p);

		p->seen=true;
		
		
		if(p->parent != NULL) {
			//sauvegarde du rang du noeud dans l'arbre de recouvrement minimum
			p->rank=p->parent->rank + 1;
			//maj des degrés
			p->parent->degre++;
			p->degre++;
			
			p->parent->fils.push_back(p);
		} 
		//sinon, c'est que c'est la racine du MST
		else if(added.size()==1){
			p->rank=0; //les ranks sont initialisé à -1
			
		}
		else {
			cout<<"Error in Prim : "<<endl<<"found a node with no parent"<<endl;
			exit(1);
		}
		
		//stop when n nodes have been added to the mst.
		if(added.size() >= graphe->getSize() ) {
			break;
		}

		//Parcourir les voisins
		 for(     nodeEdges::const_iterator it(p->nodeRef->getEdges()->begin()),
	     itstop(p->nodeRef->getEdges()->end()); it!=itstop; ++it) {
	    	
	    	int sourceNodeIndex = (*it)->getSource()->getIndex();
	    	int destNodeIndex = (*it)->getDest()->getIndex();
	    	int vref;

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
	    	Pointeur * voisinPointeur = pointeurs2[vref];
	    	
	    	double weight = (*it)->getWeight();
	    	
	    	weight+=getPenalite(sourceNodeIndex,destNodeIndex);

	    	if( !voisinPointeur->seen) {
				
				cutSize++;
				
	    		//Si le voisin est plus proche
	    		if(voisinPointeur->weight >= weight) {
	    			//on change sa distance
		    		voisinPointeur->parent = p;
		    		voisinPointeur->weight=weight;
		    		voisinPointeur->realWeight=weight;
		    		voisinPointeur->edge=*it;
		    	}
		    	
		    	//Si on est sur un arc forcé, on le met au minimum du heap
		    	if((*it)->isForced()) {
		    		voisinPointeur->parent = p;
		    		voisinPointeur->weight=-std::numeric_limits<double>::infinity();
		    		voisinPointeur->realWeight=weight;
		    		voisinPointeur->edge=*it;
		    	}
		    	//update the heap
		    	//if the neighboor is not in the heap, put it in
				if(voisinPointeur->tabIndex== -1) {
					heap.push_heap(voisinPointeur);
				}
				else { //else decrease key
					heap.monterDsHeap(voisinPointeur);
				}
			}
			else {
				cutSize--;
			}
    	}
		if(debug) {
			cout<<added.size()<<" "<<heap.getTab()->size()<<"  "<<cutSize<<endl;
		}
		//if found a cut of size <2, then the graph cannot contain a tour
		if(cutSize <=1) {
			if(debug) {
				cout<<"Found a cut of size "<<cutSize<<endl;
			}
			*containsTour = false;
			break;
		}
	
		//if there is only one node in the heap, and more than one node left to add
		//the graph cannot contains a tour
		if(added.size() <= graphe->getSize()-2 && heap.getTab()->size() == 1) {
			if(debug) {
				cout<<"There is only one node in the heap and  "<<added.size()<<" nodes in the MST"<<endl;
			}
			*containsTour = false;
			break;
			
		}
		
		
		
		//if found a cut of size 2, then edges in this cut can be forced
		if(cutSize ==2) {
			
			//if there is exactly 2 nodes in the heap
			//edges connecting theses two nodes are the two edges of the cut
			if(heap.getTab()->size() == 2) {
				//TODO : check if these edges has not already been added to the list of edges to force
				for(std::vector< Pointeur*>::iterator it = heap.getTab()->begin(); it != heap.getTab()->end(); ++it) {
					if((*it)->edge->isForced() == false) {
						toForce->push_back((*it)->edge);
					}
				}
			}
			//if there is only one node in the heap,
			//then this node have only two edges which form the cut of size 2
			else if(heap.getTab()->size() == 1) {
				Node * node = (*(heap.getTab()->begin()))->nodeRef;
				
				if(node->getEdges()->size() !=2) {
					cout<<"Error while forcing edges using cuts."<<endl;
					cout<<"Should find a node with only two edges, whereas has "<<node->getEdges()->size()<<endl;
					exit(1);
				}
				for(nodeEdges::iterator it = node->getEdges()->begin(); it != node->getEdges()->end(); ++it) {
					if((*it)->isForced() == false) {
						toForce->push_back((*it));
					}
				}
				
			}
			else {
				cout<<"Error while forcing edges using cuts."<<endl;
				cout<<"Found a cut of size 2 with "<<heap.getTab()->size()<<" nodes in the heap."<<endl;
				exit(1);
			}
		}
		
	
	}

	//free memory of pointeurs
	for(vector<Pointeur*>::iterator it = pointeurs2.begin(); it!= pointeurs2.end(); ++it) {
		delete (*it);
	}
	
	return ;
	
}



