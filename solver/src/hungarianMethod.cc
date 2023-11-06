#include <hungarianMethod.h>

extern bool verbose;

hungarianMethod::hungarianMethod(SymTSP* stsp) :
	isFromATSP(stsp->isFromATSP),
	matchingSize(0),
	matchingCost(0.0)
{
	bigraphe = new BiGraphe(stsp);
	
}

hungarianMethod::~hungarianMethod() {
	delete bigraphe;
};




void hungarianMethod::initPotential(bool * canContainTour) {
	
	//for each right node, find the smallest edge connected to it
	//the weight of this edge will be the initial potential for the node
	for(list<BiNode*>::iterator it = rightNodes_notInTree.begin(); it != rightNodes_notInTree.end(); ++it) {
		double min = std::numeric_limits<double>::infinity();
		for(bigraphEdges::iterator e_it = (*it)->getEdges()->begin(); e_it != (*it)->getEdges()->end(); ++e_it) {
			if((*e_it)->getWeight() < min) {
				min = (*e_it)->getWeight();
			}
		}
		
		if((*it)->getEdges()->size() == 0) {
			if(verbose) cout<<"Error while initializing potentials for hungarian method"<<endl;
			if(verbose) cout<<"Found a node with no edge. "<<(*it)->print()<<endl;
			(*canContainTour) = false;
			return;
			//exit(1);
		}
		
		(*it)->setPotential(min);
	}
	
	
	
	
}


double hungarianMethod::getBound() {
		return matchingCost;
}

void hungarianMethod::computeBound(bool * canContainTour) {
	if(verbose) cout<<endl<<endl<<"Starting hungarian method"<<endl;
	binodeEdges * toRemove = new binodeEdges();
	if(verbose) cout<<"graphe size : "<<bigraphe->getEdges()->size()<<endl;
	init(toRemove);
	
	initPotential(canContainTour);
	if( (*canContainTour) == false) {
		//put back in bigraphe edges removed while initialization
		for(binodeEdges::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
			bigraphe->addEdge( (*it) );
		}
		if(verbose) cout<<"graph cannot contain matching"<<endl;
		delete toRemove;
		return;
	}
	
	if(verbose) cout<<"graphe size : "<<bigraphe->getEdges()->size()<<endl;
	if(verbose) cout<<"initial matching : "<<matchingSize<<" "<<matchingCost<<endl;
	for(list<BiNode*>::iterator it = leftNodes_toMatch.begin(); it != leftNodes_toMatch.end(); ++it) {
		BiNode * node_endOfPath = NULL;
		
		if( (*it)->getMatchingEdge() != NULL) {
			if(verbose) cout<<"Error. Found a right node already matched when it should not be."<<endl;
			exit(1);
		}

		node_endOfPath =findAlternatePath((*it), canContainTour );
		
		if( (*canContainTour) == false) {
			if(verbose) cout<<"graph cannot contain matching"<<endl;
			//put back in bigraphe edges removed while initialization
			for(binodeEdges::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
				bigraphe->addEdge( (*it) );
			}
			delete toRemove;
			return;
		}
		
		if(node_endOfPath == NULL) {
			cout<<"Error. Cannot find an alternate path"<<endl;
			exit(1);
		}

		switchPath( (*it), node_endOfPath );
		
		matchingSize++;
		
	}
	
	
	for(list<BiNode*>::iterator it = leftNodes_toMatch.begin(); it != leftNodes_toMatch.end(); ++it) {
		if((*it)->getMatchingEdge() == NULL) {
			cout<<"Error. Matching not complete."<<endl;
			exit(1);
		}
		matchingCost+=(*it)->getMatchingEdge()->getWeight();
		
	}
	if(verbose) cout<<"Matching cost "<<matchingCost<<"  "<<matchingSize<<endl;
		
	
	//put back in bigraphe edges removed while initialization
	for(binodeEdges::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
		bigraphe->addEdge( (*it) );
	}
	
	
	delete toRemove;
}


BiNode * hungarianMethod::findAlternatePath(BiNode * startNode , bool * canContainsTour) {
	*canContainsTour = true;
	
	if(startNode->isInLeftSet() == false) {
		cout<<"Try to find an alternate path starting from a right node."<<endl;
		exit(1);
	}
	
	//reinitialize parents to NULL
	//slacks to infinity
	//clear slack edges
	for(bigraphNodes::iterator it = bigraphe->getRightNodes()->begin(); it != bigraphe->getRightNodes()->end(); ++it) {
		(*it)->clear();
	}
	for(bigraphNodes::iterator it = bigraphe->getLeftNodes()->begin(); it != bigraphe->getLeftNodes()->end(); ++it) {
		(*it)->clear();
	}
	
	//clear tree lists
	leaves.clear();
	leftNodes_inTree.clear();
	
	for(std::list<BiNode*>::iterator it = rightNodes_inTree.begin(); it != rightNodes_inTree.end(); ++it) {
		rightNodes_notInTree.push_back( (*it) );
		(*it)->setPosition(--(rightNodes_notInTree.end()));
	}
	rightNodes_inTree.clear();
	
	
	bool matchingFound = false;
	
	BiNode * node_endOfPath = NULL;

	//add the start node to the tree
	node_endOfPath = expandAlternateTree(startNode, &matchingFound);
	if(matchingFound) {
		return node_endOfPath;
	}
	
	
	while(matchingFound == false) {
		
		if(leaves.size()==0 ) {
			//modify potential

			
			
			//find min slack and edges corresponding to this minimum
			double min = std::numeric_limits<double>::infinity();
			std::vector<BiEdge*> minSlackEdges;
			for(std::list<BiNode*>::iterator it = rightNodes_notInTree.begin(); it != rightNodes_notInTree.end(); ++it) {
				if((*it)->getSlack() < min - PRECISION ) {
					min = (*it)->getSlack();
					minSlackEdges.assign((*it)->getMinSlackEdges()->begin(), (*it)->getMinSlackEdges()->end() );
				}
				else if( (*it)->getSlack() >= min - PRECISION &&
								(*it)->getSlack() <= min + PRECISION ) {
					minSlackEdges.insert(minSlackEdges.end(),
					(*it)->getMinSlackEdges()->begin(), (*it)->getMinSlackEdges()->end() );		
				}
				
			}
			
			for(bigraphEdges::iterator it = bigraphe->getEdges()->begin(); it != bigraphe->getEdges()->end(); ++it) {
				(*it)->getReducedCost();
			}
			
			//modify potential
			for(std::list<BiNode*>::iterator it = rightNodes_inTree.begin(); it != rightNodes_inTree.end(); ++it) {
				(*it)->setPotential((*it)->getPotential() -min);
			}
			for(std::list<BiNode*>::iterator it = leftNodes_inTree.begin(); it != leftNodes_inTree.end(); ++it) {
				(*it)->setPotential((*it)->getPotential() +min);
			}
		
			//add new leaves to the tree
			for(std::vector<BiEdge*>::iterator it = minSlackEdges.begin(); it!= minSlackEdges.end(); ++it) {
				BiNode * rightNode = (*it)->getRightNode();
				BiNode * leftNode = (*it)->getLeftNode();
				
				rightNode->setParentEdge((*it));
				if(rightNode->getMatchingEdge() == NULL) {
					matchingFound = true;
					
					return rightNode;
				}
				
				//there may be two min slack edges going to the same right node,
				//and we don't want to add it two times
				if(rightNode->isInTree()) {
					continue;
				}
				leaves.push_back(rightNode);
				
				rightNode->putInTree();
				rightNodes_notInTree.erase(rightNode->getPosition() );
				rightNodes_inTree.push_back(rightNode);
			}
			//if change of potential has not inserted a new leaf in the tree,
			//if means that a matching cannot be found
			if(leaves.size() == 0 ) {
				*canContainsTour = false;
				return node_endOfPath;
				/*
				cout<<"Change of potential has not inserted a new leaf in the tree."<<endl;
				cout<<"delta "<<min<<endl;
				cout<<"nb of min slack edges "<<minSlackEdges.size()<<endl;
				cout<<"rightNodes not in tree "<<rightNodes_notInTree.size()<<endl;
				for(std::list<BiNode*>::iterator it = rightNodes_notInTree.begin(); it != rightNodes_notInTree.end(); ++it) {
					cout<<(*it)->print();
					cout<<" ";
				}
				cout<<endl;
				cout<<"left nodes in tree "<<leftNodes_inTree.size()<<endl;
				for(std::list<BiNode*>::iterator it = leftNodes_inTree.begin(); it != leftNodes_inTree.end(); ++it) {
					cout<<(*it)->print();
					cout<<" ";
				}
				cout<<endl;
				cout<<"rightNodes in tree "<<rightNodes_inTree.size()<<endl;
				for(std::list<BiNode*>::iterator it = rightNodes_inTree.begin(); it != rightNodes_inTree.end(); ++it) {
					cout<<(*it)->print();
					cout<<" ";
				}
				cout<<endl;
				cout<<"matching size "<<matchingSize<<endl;
				printDot("uh.dot");
				bigraphe->print();
				exit(1);
				* */
			}
		}
		
		//increase the size of the tree
		while(leaves.size() >0) {
			BiNode * rightNode = leaves.front();
			leaves.pop_front();
			//all leaves of the tree are  matched 
			if(rightNode->getMatchingEdge() == NULL) {
				cout<<"Error. Found a leaf which is not matched."<<endl;
			}
			BiNode * leftNode = rightNode->getMatchingEdge()->getLeftNode();

			node_endOfPath =expandAlternateTree(leftNode, &matchingFound);
			if(matchingFound) {
				return node_endOfPath;
			}
		}
	}
	
	
	return NULL;
}

BiNode* hungarianMethod::expandAlternateTree(BiNode * node, bool * matchingFound) {

	if(node->isInLeftSet() == false) {
		cout<<"Try to expand alternate tree from a right node."<<endl;
		exit(1);
	}
	
	if(node->isInTree() ) {
		cout<<"Adding to the tree a node already in the tree."<<endl;
		exit(1);
	}
	
//	cout<<"expand tree from node "<<node->print()<<endl;
	
	leftNodes_inTree.push_back(node);
	node->putInTree();
	
	//put neigboor of node as leaves.
	for(bigraphEdges::iterator it = node->getEdges()->begin(); it != node->getEdges()->end(); ++it) {
		
		BiNode * rightNode = (*it)->getRightNode();
		
		//if edge is not tight, update slack
		if((*it)->getReducedCost() > PRECISION) {
		//	cout<<(*it)->print()<<" is not tight "<<(*it)->getReducedCost()<<endl;
		//	cout<<"right node slack "<<rightNode->getSlack()<<endl;
			if(rightNode->getSlack() > (*it)->getReducedCost() + PRECISION ) {
				rightNode->addMinSlackEdge( (*it) );	
			}
			else if ( rightNode->getSlack() <= (*it)->getReducedCost() + PRECISION && 
							rightNode->getSlack() >= (*it)->getReducedCost() - PRECISION ) {
				rightNode->addMinSlackEdge( (*it) );
			}
			continue;
		}
		
	//	cout<<(*it)->print()<<" is tight "<<(*it)->getReducedCost()<<endl;
		//don't look at right nodes already in the tree
		if(rightNode->getParentEdge() != NULL) {
			continue;
		}
		
		
		rightNode->setParentEdge( (*it));
		//stop if a free right node is found
		if(rightNode->getMatchingEdge() == NULL ) {
			*matchingFound = true;
			//(node_endOfPath) = rightNode;
			return rightNode;
		}
		//store nodes in the tree
		leaves.push_back( rightNode );
		rightNode->putInTree();
		
		rightNodes_notInTree.erase(rightNode->getPosition() );
		rightNodes_inTree.push_back(rightNode);
	}
	
	
		
	return NULL;
}

void hungarianMethod::switchPath(BiNode * startNode, BiNode * node_endOfPath) {
	
	if(node_endOfPath->isInLeftSet() ) {
		cout<<"Error. Try to switch a path ending at a left node."<<endl;
		exit(1);
	}
	
	if(node_endOfPath->getMatchingEdge() != NULL) {
		cout<<"Error. Try to switch a path ending at a matched node."<<endl;
		exit(1);
	}
	if(node_endOfPath->getParentEdge() == NULL) {
		cout<<"Error. Try to switch a path ending at a node with no parent."<<endl;
		exit(1);
	}
	
	

	
	BiNode * rightNode = node_endOfPath;
	BiEdge * parentEdge = rightNode->getParentEdge();
	BiNode *  leftNode = parentEdge->getLeftNode();
	BiEdge * old_matchingEdge = leftNode->getMatchingEdge();
	
	rightNode->setMatchingEdge(parentEdge);
	leftNode->setMatchingEdge(parentEdge);
	

	while(old_matchingEdge != NULL) {
		//cout<<rightNode->print()<<endl;
		rightNode = old_matchingEdge->getRightNode();
		parentEdge = rightNode->getParentEdge();
		
		if(parentEdge == NULL) {
			cout<<"Error while switchint a path. Found a right node with no parent"<<endl;
			exit(1);
		}
		leftNode = parentEdge->getLeftNode();
		old_matchingEdge = leftNode->getMatchingEdge();
		
		
		rightNode->setMatchingEdge(parentEdge);
		leftNode->setMatchingEdge(parentEdge);
		
	}
	
}


void hungarianMethod::init(binodeEdges * toRemove) {
	
	
	
	//init slacks and matching
	for(bigraphNodes::iterator it = bigraphe->getLeftNodes()->begin(); it != bigraphe->getLeftNodes()->end(); ++it) {
		(*it)->init();
	}
	
	for(bigraphNodes::iterator it = bigraphe->getRightNodes()->begin(); it != bigraphe->getRightNodes()->end(); ++it) {
		(*it)->init();
	}
	
	matchingCost = 0.0;
	matchingSize = 0;
	
	
	//clear lists
	leftNodes_toMatch.clear();
	rightNodes_notInTree.clear();
	leftNodes_inTree.clear();
	rightNodes_inTree.clear();
	
	//put right nodes in the list "not in the tree"
	for(bigraphNodes::iterator it = bigraphe->getRightNodes()->begin(); it != bigraphe->getRightNodes()->end(); ++it) {
		
		
		if( (*it)->getMasterNode()->getNbForcedVoisins() >2) {
			cout<<"Error while initializing hungarian method."<<endl;
			cout<<"Found a node with 3 forced edges."<<endl;
			exit(1);
		}
		//don't put nodes with forced edges
		if( (*it)->getMasterNode()->getNbForcedVoisins() ==2) {
			toRemove->insert(toRemove->end(), (*it)->getEdges()->begin(), (*it)->getEdges()->end() );
			
			//find the forced edge and add its cost to the matching
			int nbForcedEdge = 0;
			for(bigraphEdges::iterator e_it =(*it)->getEdges()->begin(); e_it != (*it)->getEdges()->end(); ++e_it) {
				if((*e_it)->getEdge()->isForced() ) {
					matchingCost+=(*e_it)->getWeight();
					matchingSize++;
					nbForcedEdge++;
					
					(*it)->setMatchingEdge( (*e_it) );
					(*e_it)->getLeftNode()->setMatchingEdge( (*e_it) );
				}
			}
			if(nbForcedEdge != 1) {
				cout<<"Error while initializing hungarian method."<<endl;
				cout<<"Found a node with two edge forced in master graph and "<<nbForcedEdge<<" forced in bigraphe."<<endl;
				exit(1);
			}
			
			(*it)->setPotential(-std::numeric_limits<double>::infinity());
			
			continue;
		}
		
		
		rightNodes_notInTree.push_back( (*it) );
		(*it)->setPosition( --(rightNodes_notInTree.end()));
	}
	
	
	//put left nodes in the list "not in the tree"
	for(bigraphNodes::iterator it = bigraphe->getLeftNodes()->begin(); it != bigraphe->getLeftNodes()->end(); ++it) {
		if( (*it)->getMasterNode()->getNbForcedVoisins() >2) {
			cout<<"Error while initializing hungarian method."<<endl;
			cout<<"Found a node with 3 forced edges."<<endl;
			exit(1);
		}
		//don't put nodes with forced edges
		if( (*it)->getMasterNode()->getNbForcedVoisins() ==2) {
			toRemove->insert(toRemove->end(), (*it)->getEdges()->begin(), (*it)->getEdges()->end() );
			(*it)->setPotential(-std::numeric_limits<double>::infinity());
			continue;
		}
		leftNodes_toMatch.push_back( (*it) );
	}
	
	
	//remove edges connected to removed nodes
	for(binodeEdges::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
		bigraphe->removeEdge( (*it) );
	}
	
	
	
}

void hungarianMethod::filter(double upperBound, double lowerBound, std::vector<Edge*>* toRemove) {
	int nb = 0;
	for(bigraphEdges::iterator it = bigraphe->getEdges()->begin(); it != bigraphe->getEdges()->end(); ++it) {
		if( (*it)->getEdge()->isForced() ) 
			continue;
		if( lowerBound + matchingCost + (*it)->getReducedCost() > upperBound + PRECISION) {
			toRemove->push_back( (*it)->getEdge() );
			nb++;
		}
	}
	if(verbose) cout<<"Filter AP : "<<nb<<" edges filtered."<<endl;
	
}

void hungarianMethod::removeEdges(std::vector<Edge*>* toRemove) {
	
	if( isFromATSP == false) {
		return;
	}
	
	for(std::vector<Edge*>::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
		bigraphe->removeEdge( (*it)->getBiEdge() );
	}
	if(verbose) cout<<"nb of edges in bigraphe "<<bigraphe->getEdges()->size()<<endl;
}

void hungarianMethod::removeEdge(Edge* edge) {
	if( isFromATSP == false) {
		return;
	}
	bigraphe->removeEdge( edge->getBiEdge() );
	if(verbose) cout<<"nb of edges in bigraphe "<<bigraphe->getEdges()->size()<<endl;
}

void hungarianMethod::addEdges(std::vector<Edge*>* toRemove) {
	if( isFromATSP == false) {
		return;
	}
	for(std::vector<Edge*>::iterator it = toRemove->begin(); it != toRemove->end(); ++it) {
		bigraphe->addEdge( (*it)->getBiEdge() );
		
	}
	if(verbose) cout<<"nb of edges in bigraphe "<<bigraphe->getEdges()->size()<<endl;
}

void hungarianMethod::addEdge(Edge* edge) {
	if( isFromATSP == false) {
		return;
	}
	
	bigraphe->addEdge( edge->getBiEdge() );
	
}

void hungarianMethod::printDot(string nomFichier) {
	bigraphe->printDot(nomFichier);
}

void hungarianMethod::setEdgesCostsInHK() {
	std::vector<double> penalites(2* bigraphe->getLeftNodes()->size(),0.0);
	
	for(bigraphEdges::iterator it = bigraphe->getEdges()->begin(); it != bigraphe->getEdges()->end(); ++it) {
		if((*it)->getEdge()->isForced() ) {
			(*it)->getEdge()->setWeight(0.0);
		}
		else {
			(*it)->getEdge()->setWeight( (*it)->getReducedCost() );
		}
	}

}

void hungarianMethod::getEdgesToBranchOn(std::vector<Edge*>* edgesToBranchOn) {
	edgesToBranchOn->clear();
	//get edges in the matching which are not forced
	for(bigraphNodes::iterator it = bigraphe->getRightNodes()->begin(); it != bigraphe->getRightNodes()->end(); ++it) {
		if((*it)->getMatchingEdge() == NULL) {
			cout<<"Error while getting matching edges: found a node not matched."<<endl;
			exit(1);
		}
		Edge * edge = (*it)->getMatchingEdge()->getEdge();
		
		//don't branch on forced edges
		if( edge->isForced() == true) {
			continue;
		}
		
		//compute replacement costs
		BiNode * leftNode = (*it)->getMatchingEdge()->getLeftNode();
		BiNode * rightNode = (*it)->getMatchingEdge()->getRightNode();
		double min = std::numeric_limits<double>::infinity();
		for(bigraphEdges::iterator e_it = rightNode->getEdges()->begin(); e_it != rightNode->getEdges()->end(); ++e_it) {
			if( (*e_it)->getLeftNode()->getIndex() == leftNode->getIndex() ){
				continue;
			}
			if((*e_it)->getReducedCost() < min) {
				min = (*e_it)->getWeight();
			}
		}
		for(bigraphEdges::iterator e_it = leftNode->getEdges()->begin(); e_it != leftNode->getEdges()->end(); ++e_it) {
			if( (*e_it)->getRightNode()->getIndex() == rightNode->getIndex() ){
				continue;
			}
			if((*e_it)->getReducedCost() < min) {
				min = (*e_it)->getWeight();
			}
		}
		
		edge->setReplacementCost(min);
		edgesToBranchOn->push_back(edge);
	}
}

bool hungarianMethod::isTour() {
	BiNode * startNode = bigraphe->getLeftNodes()->front();
	int rightNodeIndex = -1; // startNode->getMatchingEdge()->getRightNode()->getIndex();
	BiNode * leftNode = startNode; //bigraphe->getLeftNodes().at(rightNodeIndex);
	int size = 0;
	while( size < bigraphe->getLeftNodes()->size() && rightNodeIndex != startNode->getIndex() ) {
		rightNodeIndex = leftNode->getMatchingEdge()->getRightNode()->getIndex();
		leftNode = bigraphe->getLeftNodes()->at(rightNodeIndex);
		size++;
	}
	if(size == bigraphe->getLeftNodes()->size() ) {
		return true;
	}
	else {
		return false;
	}
}
