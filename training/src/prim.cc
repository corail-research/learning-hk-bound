#ifndef PRIM_CC
#define PRIM_CC
#include <stsp.h>
#include <heap.h>

SymTSP * PrimMieux() {

	bool debug = false;
	weightedNode* added[getSize()]; //Ens des noeuds de l arbre
	std::vector<weightedNode*> toAdd;
	toAdd.reserve(getSize());
	weightedNode* wnodes[maxSize]; // tableau pour acceder au wnodes
	
	//initialisation
	for(std::map<int,Node>::iterator it = nodes.begin(); it!= nodes.end(); ++it) {
		

		weightedNode* wNode = new weightedNode(&(it->second),std::numeric_limits<double>::infinity());
		wnodes[it->first]=wNode;
		toAdd.push_back(wNode);
		
	}
	
	if(debug) {
	  	cout << "myvector contains:";
		for (unsigned i=0; i<toAdd.size() ; i++)
			cout << " " << toAdd[i]->nodeRef->getIndex();

		cout << endl;
	}
	
//	added[0] =  weightedNode( &(nodes[0]), 0.0);
	int lastAddIndex=-1;
	while(toAdd.size() > 0) {
		lastAddIndex++;
	    //Recuperer le noeud le plus proche
	    
	    weightedNode* lastAdd = toAdd.front();
	    //On l'enleve des noeuds à visiter
	   pop_heap (toAdd.begin(), toAdd.end(), wnodeCompare());
	   toAdd.pop_back();
	    
	    //On le marque commme vu
	    added[lastAddIndex]=lastAdd;
		lastAdd->seen=true;
	    
	    
	    //Mettre a jour les voisins
	    Node* lastAddNode = lastAdd->nodeRef;
	  	  if(debug) {
	   		 cout <<"noeud "<< lastAddNode->getIndex() <<" added in position "<< lastAddIndex
	   		 <<" weight "<< lastAdd->weight<<"  parent "<<lastAdd->parent<<endl;
	    }
	    bool voisModified=false;
	    int countVoisMod=0;
	    int countVoisNotSeen=0;
	    int count=0;

	    for(     std::set<Edge>::const_iterator it(lastAddNode->getEdges()->begin()),
	     itstop(lastAddNode->getEdges()->end()); it!=itstop; ++it) {
	    	count++;
	    	int sourceNodeIndex = it->getSource()->getIndex();
	    	int destNodeIndex = it->getDest()->getIndex();
	    	int vref;
	    	//it->print();
	    	if(sourceNodeIndex==lastAddNode->getIndex()) {
	    		vref=destNodeIndex;
	    	}
	    	else if (destNodeIndex==lastAddNode->getIndex()) {
	    		vref=sourceNodeIndex;
	    	}
	    	else {
	    		cout << "Error in edges of node "<<lastAddNode->getIndex() <<endl
	    			<<"Found an edge from "<<sourceNodeIndex <<" to "<<destNodeIndex<<endl;
    			exit(1);
	    	}
	    	
	    	//vérifier s'il existe...
	    	//ATTENTION se fait en log(size) si wnodes est une map
	    	weightedNode* voisWnode = wnodes[vref];
	    	

	    	if(!voisWnode->seen) {
	    		countVoisNotSeen++;
	    		double edgeWeight = it->getWeight();
	    		if( voisWnode->weight > edgeWeight) {
	    			countVoisMod++;
	    			//On change la valeur de distance si il faut
	    			voisWnode->parent=lastAddNode->getIndex();
	    			voisWnode->weight=edgeWeight;
	    			voisModified=true;
	    			//push_heap(toAdd.begin(),positionDevoisWnodeDanstoAdd, wnodeCompare()):
	    			if(debug) {
						cout << "node " <<voisWnode->nodeRef->getIndex()<< " newparent "<<voisWnode->parent
						<<" weight "<<edgeWeight<<endl;
	    			}			
	    		}
	    	}

	    	
	    }

	    //mettre à jour le heap
	    if(voisModified) {
		   	 make_heap(toAdd.begin(), toAdd.end(), wnodeCompare());
		   	 
			if(debug) {
			  	cout << "after maj voisins:";
				for (unsigned i=0; i<toAdd.size() ; i++)
					cout << " " << toAdd[i]->nodeRef->getIndex();

				cout << endl;
			}
	    }
	    
	}
	
	SymTSP * arbre = new SymTSP(getSize());

	for(int i=1; i< getSize(); i++) {
		
		int source = (*(added+i))->nodeRef->getIndex();
		int dest = (*(added+i))->parent;
		double weight =(*(added+i))->weight;
		if(debug) {
			cout << "node "<<source<<" has parent "<<dest<<endl;
		}
		arbre->addEdge(source, dest, weight);
	}
	

	
	return arbre;
}

#endif
