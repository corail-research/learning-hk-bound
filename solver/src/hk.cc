#include <hk.h>
#include <Python.h>
#include <float.h>
#include <chrono>

bool debug=false;
extern bool heavyFiltering;
extern bool verbose;

bool callPythonThetas = true;
int cHKBnB = 0;
int HKBnB_limit = 10;

HeldKarp::HeldKarp(SymTSP * g) :
	graphe(g), 
	computed(false),
	oneEdge1(NULL),
	oneEdge2(NULL)
	{
		//min1tree = dynamic_cast<SymTSP*>(g->getSameGraphWithoutEdges());
		int n= g->getMaxSize();
		
		dfscount=0;
		//initialisation des pointeurs servant dans Prim.
		pointeurs=vector<Pointeur*>(n);

		int i=0;
		
		//reserve la mémoire pour les pénalités.
		penalite=vector<double>(n,0.0);
		bestPenalites=vector<double>(n,0.0);
		
		degrees =vector<int>(n,0);
		bestDegrees =vector<int>(n,0);
		
		
		nodeExist = vector<bool>(n,false);
		
//		oneTreeEdges.reserve(n);
		

		
		parents=vector<parentMaxEdge>(n);
		
		for(graphNodes::iterator it = graphe->getNodes()->begin(); it!= graphe->getNodes()->end(); ++it) {
		
			//Attention à bien initialiser les pointeurs avec leur position dans le heap
			//Pointeur * p = new Pointeur(&(it->second), i, std::numeric_limits<double>::infinity());
			//i++;

				pointeurs.at((*it)->getIndex())=new Pointeur((*it), i, std::numeric_limits<double>::infinity());
				i++;


			nodeExist.at((*it)->getIndex())=true;
		}
		
	
		//cout<<"pointeurs.size : "<<pointeurs.size()<<endl;
		
	}
	
	
	
double HeldKarp::getBound(){
	return HKBound;
}
	
double HeldKarp::getPenalite(int source, int dest) {
	if(source==dest) {
		cout << "cannot find penalite for edge "<<source<<" to " << dest<<endl;
		exit(1);
	}
	
	if(source <0 || source > graphe->getMaxSize() || !nodeExist.at(source)) {
		cout << "cannot find penalite for edge "<<source<<" to " << dest<<endl;
		cout << "There is no node "<<source<<endl;
		exit(1);
	}
	
	if(dest <0 || dest > graphe->getMaxSize()  ||  !nodeExist.at(dest) ){
		cout << "cannot find penalite for edge "<<source<<" to " << dest<<endl;
		cout << "There is no node "<<dest<<endl;
		exit(1);
	}
	
	/*
	int i,j;
	
	
	if(source <dest) {
		i=source;
		j=dest;
	}
	else {
		i=dest;
		j=source;
	}
	
	int maxSize = min1tree->getMaxSize();
	int tabIndex=i*maxSize-(int)((i*i+i)/2.0);
	*/
	
	
	return penalite.at(source)+penalite.at(dest);
}

void HeldKarp::dfs(Pointeur* source, parentMaxEdge pme) {


//	cout<<"dfs on "<<source->nodeRef->getIndex()<< endl
	//<<"degre : "<<source->degre<<endl
		//	<<"fils.size : "<<source->fils.size()<<endl;
	//on regarde si l'arc montant est le plus grand
	if(source->parent !=NULL) {
		if(source->weight > pme.weight) {
			pme.maxEdgeSourceIndex=source->nodeRef->getIndex();
			pme.weight=source->weight;
		}
		
		if(pme.parentIndex==pme.maxEdgeSourceIndex) {
			pme.maxEdgeSourceIndex=source->nodeRef->getIndex();
			pme.weight=source->weight;
		}
	}
	//on sauvegarde le couple parent/arc montant max
	parents[source->nodeRef->getIndex()]=pme;
	//cout<<"parent de "<<source->nodeRef->getIndex()<<" set"<<endl;
	dfscount++;
	
//Si on est sur une feuille, on arrete
	if(source->fils.size()==0) {
		return;
	}
	
//Si on est sur un noeud parent, on réinitialise le pme
	if(source->degre >2) {
	//if(source->fils.size() >1 ) {
		pme.weight=-std::numeric_limits<double>::infinity();
		pme.parentIndex=source->nodeRef->getIndex();
		pme.maxEdgeSourceIndex=source->nodeRef->getIndex();
	}
	
	bool changeNumBranche=false;
	
	if(source->nodeRef->getIndex() == index1Edge1 || source->nodeRef->getIndex() == index1Edge2 && source->degre ==3) {
		changeNumBranche=true;
	}
	
	if(source->degre==2 && source->parent!=NULL) {
		if(source->fils.size()!=1) {
			cout<<"Erreur dans dfs, sommet de degré 2 qui n'est pas la racine et qui a plus d'un fils"<<endl
			<<"sommet : "<<source->nodeRef->getIndex()<<endl
			<<"degre : "<<source->degre<<endl
			<<"fils.size : "<<source->fils.size()<<endl;
			
			for(std::vector<Pointeur*>::iterator it = source->fils.begin();
					it!=source->fils.end(); ++it) {
				cout<<(*it)->nodeRef->getIndex()<<"  ";
			}
			cout<<endl;
			exit(1);
		}
	}
	
//On parcourt les enfants
	int branche=0;
	for(std::vector<Pointeur*>::iterator it = source->fils.begin();
				it!=source->fils.end(); ++it) {
		if(source->fils.size()>1 || changeNumBranche) {
			pme.branche=branche;
		}
		dfs(*it,pme);			
		branche++;
	}
	


}


//double : poids de l'arc avec penalite
//bool = false si l'arc est dans le 1tree
void HeldKarp::getWeightOfRemovedEdge(Edge* edge, bool * edgeIsInOneTree, double * weightOfRemovedEdge) {


	int source =edge->getSource()->getIndex();
	
	int dest = edge->getDest()->getIndex();


	//!!!!! en log(m) !!!
	/*
	if(!graphe->containsEdge(source,dest)) {
		cout <<"There is no edge from "<<source<<" to "<<dest<<". Cannot find the cycle created."<<endl;
		exit(1);
	}
	*/
	if(source==index1node) {
		
		//if edge is oneEdge1 or oneEdge2, return
		if(dest==index1Edge1) {
			*weightOfRemovedEdge  = oneEdge1Weight;
			*edgeIsInOneTree      = true;
			return;
		}
		if(dest==index1Edge2) {
			*weightOfRemovedEdge  = oneEdge2Weight;
			*edgeIsInOneTree      = true;
			return;
		}
		//here edge is  neither oneEdge1 nor oneEdge2
		
		//edge can be replaced by oneEdge2
		if(! oneEdge2->isForced() ) {
			*weightOfRemovedEdge = oneEdge2->getWeight()+getPenalite(index1node, index1Edge2);
			*edgeIsInOneTree     = false;
			return;
		}
		//edge cannot be replaced, so its reduced cost is infinity
		//thus the weight of the removed edge if -infinity
		else if(oneEdge1->isForced()) { 
			cout<<"les 2 arcs du 1node sont forcé, et il y a des arcs connecté au 1node non éliminés"<<endl;
			*weightOfRemovedEdge  = -std::numeric_limits<double>::infinity();
			*edgeIsInOneTree      = false;
			return;
		}
		//edge can be replaced by oneEdge1
		else {
			cout<<"1edge 2 est forcé et pas 1edge1..."<<endl;
			*weightOfRemovedEdge =oneEdge1->getWeight()+getPenalite(index1node, index1Edge1);
			*edgeIsInOneTree      = false;
			return;
		}	
	}

	if(dest==index1node) {
		if(source==index1Edge1) {
			*weightOfRemovedEdge  = oneEdge1Weight;
			*edgeIsInOneTree      = true;
			return;
		}
		if(source==index1Edge2) {
			*weightOfRemovedEdge  = oneEdge2Weight;
			*edgeIsInOneTree      = true;
			return;
		}
		//here edge is  neither oneEdge1 nor oneEdge2
		
		//edge can be replaced by oneEdge2
		if(!oneEdge2->isForced()) {
			*weightOfRemovedEdge = oneEdge2->getWeight()+getPenalite(index1node, index1Edge2);
			*edgeIsInOneTree     = false;
			return;
		}
		//edge cannot be replaced, so its reduced cost is infinity
		//thus the weight of the removed edge if -infinity
		else if(oneEdge1->isForced()) { 
			cout<<"les 2 arcs du 1node sont forcé, et il y a des arcs connecté au 1node non éliminés"<<endl;
			*weightOfRemovedEdge  = -std::numeric_limits<double>::infinity();
			*edgeIsInOneTree      = false;
			return;
		}
		//edge can be replaced by oneEdge1
		else {
			cout<<"1edge 2 est forcé et pas 1edge1..."<<endl;
			*weightOfRemovedEdge =oneEdge1->getWeight()+getPenalite(index1node, index1Edge1);
			*edgeIsInOneTree      = false;
			return ;	
		}	
		
	}
	
	
	//if edge is not connected to the one node, 
	//find the max edge on the cycle created by edge in the MST
	getWeightOfMaxEdgeOnCycleCreatedBy(edge, edgeIsInOneTree, weightOfRemovedEdge);

	return ;

}

void HeldKarp::getWeightOfMaxEdgeOnCycleCreatedBy(Edge* edge,  bool * edgeIsInOneTree, double * weightOfRemovedEdge) {


	int source =edge->getSource()->getIndex();
	
	int dest = edge->getDest()->getIndex();
	
	double outEdgeWeight= edge->getWeight()+getPenalite(source, dest);




	//On cherche l'extremité de l'arete qui est de plus petit rang dans l'arbre de recouvrement minimum.
	

	Pointeur* sourcePointeur = pointeurs.at(source);
	Pointeur* destPointeur = pointeurs.at(dest);
	
	Pointeur* i=destPointeur;
	Pointeur* j=sourcePointeur;
	if(i->rank < j->rank) {
		i=sourcePointeur;
		j=destPointeur;
	}

	//if i's parent is j, then edge is in the MST
	if(i->parent->nodeRef->getIndex() == j->nodeRef->getIndex() ) {

		*weightOfRemovedEdge = i->weight;
		*edgeIsInOneTree     = true;
		return;

	}

	pair<int,int> result(i->nodeRef->getIndex(), i->parent->nodeRef->getIndex() );
	
	//attention, getEdge se fait en log(m)
//	double maxWeight=graphe->getEdge(result.first, result.second)->getWeight();
	
	double maxWeight = i->weight;
	
	int nbParentCroise=0;
	
	std::vector<int> cheminI, cheminJ;
	
	//On remonte i dans l'arbre tant que rang[i] > rang[j]
	while(i->rank > j->rank ) {
	
		if(i->fils.size()>1) {
			nbParentCroise++;
			
		}
		
		pair<int,int> arcCourant(i->nodeRef->getIndex(), i->parent->nodeRef->getIndex() );
		//double weight = graphe->getEdge(arcCourant.first, arcCourant.second)->getWeight();
		double weight = i->weight;
		
		if(weight > maxWeight) {
			maxWeight=weight;
			result=arcCourant;
		}
		
		//Retenir l'arc testé dans le support de i
		i->addEdgeInSupport(edge ,outEdgeWeight );
		
		i= i->parent;
		cheminI.push_back(i->nodeRef->getIndex());
		

		
	}
	//normalement ici rank[i]=rank[j]
	if(i->rank != j->rank ) {
		cout<<"pb de rang"<<endl;
		exit(1);
	}
	
	//On remonte dans les 2 branches jusqu'à trouver un ancetre commun
	
	while( i->nodeRef->getIndex() != j ->nodeRef->getIndex() ) {
	
		//Ne devrait jamais arriver...
		if( i->parent ==NULL || j->parent == NULL ) {
			cout <<"Error while computing cycle created by edge "<<source<<" to "<<dest<<". Could not find a cycle."<<endl;
			exit(1);
		}
		
		pair<int,int> arcCourant(i->nodeRef->getIndex(), i->parent->nodeRef->getIndex() );
		//double weight = graphe->getEdge(arcCourant.first, arcCourant.second)->getWeight();
		double weight = i->weight;
		
		if(weight > maxWeight) {
			maxWeight=weight;
			result=arcCourant;
		}
		if(i->fils.size()>1) {
			nbParentCroise++;
			
			
		}
		//Retenir l'arc testé dans le support de i
		i->addEdgeInSupport(edge ,outEdgeWeight );
		
		i= i->parent;
		cheminI.push_back(i->nodeRef->getIndex());
		
		arcCourant= pair<int,int>(j->nodeRef->getIndex(), j->parent->nodeRef->getIndex() );
		//weight = graphe->getEdge(arcCourant.first, arcCourant.second)->getWeight();
		weight = j->weight;
		
		if(weight > maxWeight) {
			maxWeight=weight;
			result=arcCourant;
		}
		if(j->fils.size()>1) {
			nbParentCroise++;
			
		}
		//Retenir l'arc testé dans le support de j
		j->addEdgeInSupport(edge ,outEdgeWeight );
		
		j= j->parent;
		cheminJ.push_back(j->nodeRef->getIndex());
	
	}
	//cout<<result.first<<" to "<<result.second<<"     "<<maxWeight<< endl;
	/*
	cout<<"chemin i : ";
	for(std::vector<int>::iterator it = cheminI.begin(); it!= cheminI.end(); ++it) {
		cout<<*it<<" ";
	}
	cout<<endl<<"chemin j : ";
	for(std::vector<int>::iterator it = cheminJ.begin(); it!= cheminJ.end(); ++it) {
		cout<<*it<<" ";
	}
	cout<<endl;
	cout<<"nb Parent croisé : "<<nbParentCroise<<endl;
	*/
	//cout <<"result \t"<<result.first<<" rank "<< pointeurs.at(result.first)->rank<<"        parent "<<i->nodeRef->getIndex()<<" rank "<<i->rank<<endl;
	
	//MaxEdge = result.first;
	
	*weightOfRemovedEdge = maxWeight;
	*edgeIsInOneTree     = false;
	
	return ;

}

void HeldKarp::filter_forcedDegrees( std::vector<Edge*>* toRemove) {
	
	for(graphNodes::iterator it = graphe->getNodes()->begin(); it != graphe->getNodes()->end(); ++it) {
		if( (*it)->getNbForcedVoisins() == 2) {
			int nbForcedEdge = 0;
			for( graphEdges::iterator e_it = (*it)->getEdges()->begin(); e_it != (*it)->getEdges()->end(); ++e_it) {
				if( (*e_it)->isForced() ) {
					nbForcedEdge++;
				}
				else {
					toRemove->push_back( (*e_it) );
				}
			}
			if(nbForcedEdge != 2) {
					cout<<"Error while filtering edges based on forced degree."<<endl;
					cout<<"Found a node which should have 2 forced edges, but have "<<nbForcedEdge<<endl;
					exit(1);
			}
		}
		else if( (*it)->getNbForcedVoisins() > 2) {
			cout<<"Error. Found a node with more than 2 forced edges while filtering based on forced degree."<<endl;
			cout<<"Should have been detected while forcing edges."<<endl;
			exit(1);
		}
		
	}
	
}

void HeldKarp::filter_reducedCosts(double upperBound, std::vector<Edge*>* toRemove) {

	filterOK=true;

	if(verbose) cout<<"start filter"<<endl;
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end(); ++it) {
		(*it)->supportCycleEdges.clear();
		(*it)->minEdgeInSupport=NULL;
		(*it)->minEdgeWeight=std::numeric_limits<double>::infinity();
	}


	int seen=0;
	
	for(graphEdges::iterator it = graphe->getEdges()->begin(); it!= graphe->getEdges()->end(); ++it) {
		
		//si l'arc est forcé, on peut pas l'enlever...
		if(  (*it)->isForced()) {
			(*it)->coutMarg=0.0;
			continue;
		}
		
		//Si un sommet de l'arc a deux arcs forcés , on peut enlever cet arc 
		//(l'arc courant ne peut pas être un arc forcé ici)
		/*
		if((*it)->getSource()->getNbForcedVoisins() >=2 ) {
			if(debug) {
				cout<<"Edge from "<<(*it)->getSource()->getIndex()<<" to "<<(*it)->getDest()->getIndex()<<" can be removed because node "<<(*it)->getSource()->getIndex()<<" has "<<(*it)->getSource()->getNbForcedVoisins()<<" forced edges"<<endl;
			}
			
			toRemove->push_back(*it);
			continue;
		}
		if((*it)->getDest()->getNbForcedVoisins() >=2 ) {
			if(debug) {
				cout<<"Edge from "<<(*it)->getSource()->getIndex()<<" to "<<(*it)->getDest()->getIndex()<<" can be removed because node "<<(*it)->getDest()->getIndex()<<" has "<<(*it)->getDest()->getNbForcedVoisins()<<" forced edges"<<endl;
			}
			
			toRemove->push_back(*it);
			continue;
		}
		*/
		
		int sourceIndex=(*it)->getSource()->getIndex();
		int destIndex= (*it)->getDest()->getIndex();
		
		
		bool  edgeIsInOneTree;
		double  weightOfRemovedEdge;
		getWeightOfRemovedEdge((*it), &edgeIsInOneTree, &weightOfRemovedEdge);
		
	
		//si result.second=false, c'est que l'arc est dans le 1tree
		if(edgeIsInOneTree) {
			(*it)->coutMarg=0.0;
			continue;
		}
		
		double augmentation = (*it)->getWeight()+getPenalite(sourceIndex,destIndex)
		 - weightOfRemovedEdge;
		 
		 (*it)->coutMarg=augmentation;
		 	
		 	
		int source = (*it)->getSource()->getIndex();
		int  dest  = (*it)->getDest()->getIndex();
		
		
		if(augmentation < - 2*1e-11){
	//		cout << -DBL_EPSILON<<endl;
			cout<<"Erreur, augmentation négative lors du filtrage"<<endl;
			cout <<"Edge from "<<sourceIndex<<" to "<<destIndex<<endl;
			cout<< "poids à rajouter : "<<(*it)->getWeight()+getPenalite(sourceIndex,destIndex)<<
			"      poids à enlever : "<<weightOfRemovedEdge<<endl<<augmentation<<"  "<<endl;
			//cout<<oneEdge1Weight<<"   "<<oneEdge2Weight<<endl;
			printDot("erreur.dot");
			exit(1);
		}
	
		if(HKBound+augmentation - upperBound > PRECISION ){
			if(debug) {
				cout <<"Edge from "<<sourceIndex<<" to "<<destIndex<<
				" can be removed because HKBound+augmentation >upperBound "<<endl;
			}
			toRemove->push_back((*it));

		}
		
	}
	
//	cout<<toRemove->size()<<" edges can be removed by filtering among "<<graphe->getEdges()->size()<<endl;
	//cout<<seen<<endl;
	return;
}


void HeldKarp::PrimMieux(bool * canContainTour) {

	bool debug = false;
	std::vector<Pointeur*> added; //Ens des noeuds de l arbre
	added.reserve(graphe->getSize()-1); //ya un -1 car on prend pas le 1node
	
	std::vector<Pointeur*> toAdd; //tableau servant pour le heap
	toAdd.reserve(graphe->getSize()); //On connait deja la taille max du vecteur : on reserve la mémoire
	
	//fabrication du heap
	Heap heap(toAdd);
	
	
	//réinitialisation des pointeurs
	int k=0;
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end(); ++it) {
		(*it)->clear();
		if((*it)->nodeRef->getIndex()==index1node) {
			oneNode=(*it);
		}
	}
	
	//initialize heap with a node which is not the one node
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end(); ++it) {
		if((*it)->nodeRef->getIndex()!=index1node) {
			heap.push_heap(*it);
			
			break;
		}
	}
	
	if(debug) {
	  	cout << "heap contains: ";
		heap.print();
		cout << endl;
	}
	

	while( heap.getTab()->size() >0) {
		
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
	    	
			
			//on ne s'occupe pas des arcs incidents au 1node
	    	if(vref== index1node) {
	    		continue;
	    	}

	    	//Attention : il faudrai vérifier que ce voisin est bien un 
	    	//element du graphe, mais prendrai log(nodes.size()) opérations...
	    	Pointeur * voisinPointeur = pointeurs[vref];
	    	
	    	double weight = (*it)->getWeight();
	    	
	    	weight+=getPenalite(sourceNodeIndex,destNodeIndex);

	    	if( !voisinPointeur->seen) {
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
    	}	
	
	}
	
	//if some nodes has not been visited, then the graph is not connected
	if(added.size() < graphe->getSize()-1) {
		*canContainTour = false;
		return;
	}
	*canContainTour = true;
	//fabrication de l'arbre résultant.

	
	double totalWeight =0.0;
	
	mstRoot=added.at(0);

	//compute MST real weight
	for(int i=1; i< added.size(); i++) { 
		
		int source = added.at(i)->nodeRef->getIndex();
		int dest = added.at(i)->parent->nodeRef->getIndex();
		double weight =added.at(i)->realWeight;
		totalWeight+= weight;

		if(debug) {
			cout << "node "<<source<<" has parent "<<dest<<endl;
		}

	}
	

	
	
	//ajouter les 2 arcs minimums incidents au 1-node.
	

	oneEdge1Weight=std::numeric_limits<double>::infinity();
	oneEdge2Weight=std::numeric_limits<double>::infinity();
	oneEdge3Weight=std::numeric_limits<double>::infinity();
	for(nodeEdges::const_iterator it = oneNode->nodeRef->getEdges()->begin(); 
							it !=oneNode->nodeRef->getEdges()->end(); ++it) {
							
							
			int sourceNodeIndex = (*it)->getSource()->getIndex();
	    	int destNodeIndex = (*it)->getDest()->getIndex();
	    	
	    	int otherNodeIndex=0;
	    	if(sourceNodeIndex != index1node) {
	    		otherNodeIndex=sourceNodeIndex;
	    	}
	    	else if (destNodeIndex!=index1node) {
	    		otherNodeIndex=destNodeIndex;
	    	}
	    	else {
	    		cout<<"Error while computing oneEdges"<<endl;
	    		exit(1);
	    	}
	    	
			
			double weight = (*it)->getWeight();
			weight+=getPenalite(index1node, otherNodeIndex);
			
			if((*it)->isForced()) {
				weight=-std::numeric_limits<double>::infinity();
			}
			

			if(weight < oneEdge1Weight) {
				oneEdge3Weight=oneEdge2Weight;
				oneEdge2Weight=oneEdge1Weight;
				oneEdge2=oneEdge1;
				index1Edge2=index1Edge1;
				
			
				oneEdge1Weight = weight;
				oneEdge1=*it;//pair<int,int>(index1node, otherNodeIndex);
				index1Edge1=otherNodeIndex;


			}
			else if(weight < oneEdge2Weight) {

			
				oneEdge3Weight=oneEdge2Weight;
				oneEdge2Weight = weight;
				oneEdge2=*it;  //pair<int,int>(index1node, otherNodeIndex);
				index1Edge2=otherNodeIndex;
				
			}
			else if (weight < oneEdge3Weight) {
				oneEdge3Weight = weight;
			}
			
							
	}
	
	//le degré du 1node est 2
	oneNode->degre=2;
	
	pointeurs.at(index1Edge1)->degre++;
	pointeurs.at(index1Edge2)->degre++;

	double oneEdgesCost = oneEdge1->getWeight()+oneEdge2->getWeight()
				+getPenalite(index1node, index1Edge1) + getPenalite(index1node, index1Edge2);
				
	
	oneTreeWeight=totalWeight+oneEdgesCost;
//	cout <<"Cost of the minimum one tree "<<oneTreeWeight<<"        mst  "<<totalWeight<<endl;	
	//return arbre;

	return;
	
}


void HeldKarp::computeBound(double upperBound, bool *tourFound, bool *canStopBranching) {

	//	cout <<endl<<"Compute bound version 1"<<endl;
	/*
	if(computed) {
		cout <<"Held and Karp bound already computed : "<<HKBound<<endl;
		return;
		
	}
*/
	*canStopBranching = false;
	*tourFound=false;
	
	
	double targetBound=1.1*upperBound;

	int nbIteration=0;
	double iterFactor = 0.15;
	int maxChange=25;
	int nbIter=iterFactor*graphe->getSize();
	if(nbIter < 5) {
		nbIter=5;
	}
	double tSmall = 0.01;
	double alpha =2.0;
	double beta =0.5;
	double step=0.0;
	
	int count_no_increase=0;
	
	bool canContainTour;


                if(callPythonThetas && cHKBnB < HKBnB_limit && cHKBnB > 0){

		Py_Initialize();
		PyRun_SimpleString("import sys");
		PyRun_SimpleString("sys.path.append(\".\")"); 

	    // Import the Python module
		PyObject* pModule = PyImport_Import(PyUnicode_DecodeFSDefault("loadNN"));
		cHKBnB++;
	    if (pModule == nullptr) {
		cout << "pModule == nullptr" << endl;
		PyErr_Print();
		return;
	    }

	    // Get the reference to the Python function
	    PyObject* pFunc = PyObject_GetAttrString(pModule, "getLogits");
	    if (pFunc == nullptr || !PyCallable_Check(pFunc)) {
		PyErr_Print();
		Py_DECREF(pModule);
		return;
	    }

	    // Prepare and pass arguments
	    PyObject *pArgs, *pResult;
		pArgs = PyTuple_New(1);
		PyObject *PList = PyList_New(0);
		PyList_Append(PList, PyLong_FromLong(cHKBnB));
		std::__cxx11::list<Edge*>* my_list = graphe->getEdges(); 

		// Iterate over the list using a pointer to list iterator
		for (auto it = my_list->begin(); it != my_list->end(); ++it) {
		    Edge* edge = *it;  
			PyList_Append(PList, PyLong_FromLong(edge->getSource()->getIndex()));
			PyList_Append(PList, PyLong_FromLong(edge->getDest()->getIndex()));
			PyList_Append(PList, PyLong_FromLong(edge->getWeight()));
			PyList_Append(PList, PyLong_FromLong(edge->isForced() ? 1 : 0));
		}


		PyTuple_SetItem(pArgs, 0, PList);



			pResult = PyObject_CallObject(pFunc, pArgs);

			if (pResult == nullptr) {
				PyErr_Print();
				Py_DECREF(pArgs);
				Py_DECREF(pFunc);
				Py_DECREF(pModule);
				return;
			    }


		    if (!PyList_Check(pResult)) {
			PyErr_Print();
			Py_DECREF(pResult);
			Py_DECREF(pFunc);
			Py_DECREF(pModule);
			return;
		    }
		    Py_ssize_t listSize = PyList_Size(pResult);

		    // Process the elements of the list
		    for (Py_ssize_t i = 0; i < listSize; ++i) {
			PyObject* pItem = PyList_GetItem(pResult, i);
			if (PyFloat_Check(pItem)) {
			    double itemValue = PyFloat_AsDouble(pItem);
				if(callPythonThetas){
				penalite[i] = itemValue;}
			}
		    }


	    Py_DECREF(pArgs);
	    Py_DECREF(pResult);
	    Py_DECREF(pFunc);

	}


	if(cHKBnB == 0){cHKBnB++;}

		toRemove.clear();
	
	PrimMieux(&canContainTour);
	
	if(!canContainTour) {
		HKBound=std::numeric_limits<double>::infinity();
		if(verbose) cout<<"Graphe pas connexe"<<endl;
		*canStopBranching = true;
		return;
	}
	
	double bestBound=-std::numeric_limits<double>::infinity();
	
	
	bool endReached=false;

	
	for(int l=1; l<maxChange+1; l++) {
		for(int k=1; k< nbIter+1; k++) {
			
			nbIteration++;
			
			//calculs préalables
			double totalPenalite=0.0;
			double violation=0.0;
			bool tour=true;
			int count=0;
			//cout<<step<<endl;
			for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end() ; ++it) {
				int index = (*it)->nodeRef->getIndex();
				totalPenalite+=penalite[index];
				
				degrees[index] = (*it)->degre;
				
				violation+=(2-(*it)->degre)*(2-(*it)->degre);
				if((*it)->degre != 2) {
					count++;
					tour = false;
				}
				
			}
			HKBound = oneTreeWeight-2*totalPenalite;
			
			
		
			
			/*
			if(heavyFiltering) {
			
				std::vector<Edge*>* toRemove2 = filter(upperBound,false);
				//printDot("arbre.dot");
				graphe->removeEdges(toRemove2);
				toRemove.insert(toRemove.end(),toRemove2->begin(), toRemove2->end() );
				delete toRemove2;
			}
			*/
	
			

			//keep best bound
			if(HKBound - bestBound > PRECISION ) {
				count_no_increase=0;
				bestBound=HKBound;
				bestPenalites.assign(penalite.begin(), penalite.end() );
				bestDegrees.assign(degrees.begin(), degrees.end() );
			}
			
			//si on a un tour, on arrete
			if(tour) {
				if(verbose) cout <<"tour found !";
				if(verbose) cout<<" Weight : "<<HKBound<<endl;
				*canStopBranching = true;
				*tourFound=true;
				return ;
			}
			
			
			if(HKBound - upperBound > PRECISION ) {
				if(verbose) cout<<"LB > UB"<<endl;
				if(verbose) cout<<"hkBound : "<<HKBound<<endl;
				if(verbose) cout<<"upperbound : "<<upperBound<<endl;
				*canStopBranching = true;
				return ;
			}
			
			
			if(violation==0) {
				cout<<"Error while computing bound. Found violation==0 whereas 1tree is not a tour..."<<endl;
				exit(1);
			}
			
			//if there is a big decrease of the bound,
			//go back to the best penalties and reduce step
			/*
			if(HKBound  <  0.95*bestBound) {

				penalite.assign(bestPenalites.begin(), bestPenalites.end()) ;
				degrees.assign(bestDegrees.begin(), bestDegrees.end() );
				alpha*=beta;
				//cout<<"going back"<<endl;
			}
			*/
			//calcul du pas
			step=alpha*(targetBound-HKBound)/violation;
			//cout<<"step "<<step<<endl;

			
			if(step<0.0) {
				cout<<"!!!!!!!!! step <0"<<endl;
				cout<<"step : "<<step<<endl;
				cout<<"hkBound : "<<HKBound<<endl;
				cout<<"upperbound : "<<upperBound<<endl;
				cout<<"alpha : "<<alpha<<endl;
				cout<<"ub-hk : "<<upperBound-HKBound<<endl;
				
				exit(1);	
			}
			
			//stop if the step is too small
			if(step < tSmall) {
				endReached=true;
				break;
			}
			
			
			//mise à jour des pénalités
			for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end() ; ++it) {
			
				int index = (*it)->nodeRef->getIndex();
				//mise à jour de la pénalité
				penalite[index]+=step*(degrees[index]-2);
			}
			PrimMieux(&canContainTour);
			
			if(!canContainTour) {
				HKBound=std::numeric_limits<double>::infinity();
				*canStopBranching = true;
				return;
			}

		}
		
		
		
		if(endReached) {
			break;
		}
		
		
		if(HKBound > bestBound- 0.1) {
			count_no_increase++;
		}
		if(count_no_increase==3) {
			break;
		}
		
		alpha*=beta;
	}
	
	//On se remet sur le meilleur point trouvé
	double totalPenalite=0.0;
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end() ; ++it) {
		int index = (*it)->nodeRef->getIndex();
		//mise à jour de la pénalité
		penalite[index]=bestPenalites[index];
		totalPenalite+=penalite[index];
	}
	PrimMieux(&canContainTour);
	if(!canContainTour) {
		HKBound=std::numeric_limits<double>::infinity();
		*canStopBranching = true;
		return;
	}
	
	HKBound=oneTreeWeight-2*totalPenalite;
	
	
	if(verbose) cout<<"nb d'iterations : "<<nbIteration<<endl;
	if(verbose) cout<<"Held and Karp bound : "<<HKBound<<endl;
	return ;
}









void HeldKarp::clearPenalties() {
	penalite=std::vector<double>(graphe->getMaxSize(), 0.0);
}

bool HeldKarp::isEdgeInMST(int source, int dest) {
	if(source==index1node || dest==index1node)
		return false;
	if(source  != mstRoot->nodeRef->getIndex()) {
		if(pointeurs.at(source)->parent->nodeRef->getIndex()==dest) {
			return true;
		}
	} 
	if(dest  != mstRoot->nodeRef->getIndex()) {
		if(pointeurs.at(dest)->parent->nodeRef->getIndex()==source ) {
			return true;
		}
	}
	return false;
}


void HeldKarp::printDot(string nomFichier){
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"graph " <<"{"<<endl;
     
      
    	
		for( std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!=pointeurs.end(); ++it) {
			//Edge ed = *it;
			if((*it)->nodeRef->getIndex() == index1node) {
				continue;
			}
			if((*it)->nodeRef->getIndex() == mstRoot->nodeRef->getIndex()) {
				continue;
			}
			
			double weight=(*it)->realWeight;//-getPenalite((*it)->nodeRef->getIndex(), (*it)->parent->nodeRef->getIndex());
			int source = (*it)->nodeRef->getIndex();
			
			if((*it)->parent==NULL) {
				continue;
			}
			
			int dest = (*it)->parent->nodeRef->getIndex();
			
			
			/*
			if(deltaPSet[source] && pointeurs.at(source)->degre >2 
			//&& source !=index1Edge1 && source != index1Edge2
			) {
				weight+=deltaP[source];
				
			}
			
			if(deltaPSet[dest] && pointeurs.at(dest)->degre >2
			// && dest !=index1Edge1 && dest != index1Edge2
			 ) {
				weight+=deltaP[dest];
				
			}
			*/
			
			if((*it)->edge->isForced()) {
				weight=std::numeric_limits<double>::infinity();
			}
			fb << (*it)->nodeRef->getIndex() <<" -- " << (*it)->parent->nodeRef->getIndex()<<
				"[label = \" "<< weight<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		
oneEdge1->getWeight()+oneEdge2->getWeight()
				+getPenalite(index1node, index1Edge1) + getPenalite(index1node, index1Edge2);
				
	
		double e1=oneEdge1->getWeight()+getPenalite(index1node, index1Edge1);
		fb << index1node <<" -- " << index1Edge1<<
				"[label = \" "<< e1<<"\" ]" <<	endl;
				
		double e2=oneEdge2->getWeight()+getPenalite(index1node, index1Edge2);
		fb << index1node <<" -- " << index1Edge2<<
				"[label = \" "<< e2<<"\" ]" <<	endl;

		fb << "}"<<endl;

        // Ferme le fichier :
        fb.close();
    }


}



void HeldKarp::printDot2(string nomFichier){
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"graph " <<"{"<<endl;
     
      
    	
		for( std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!=pointeurs.end(); ++it) {
			//Edge ed = *it;
			
			fb<< (*it)->nodeRef->getIndex()<< "[ label = \" "<<(*it)->nodeRef->getIndex()<<" \\n "<<penalite.at( (*it)->nodeRef->getIndex())<<" \" ]"<<endl;
			
			if((*it)->nodeRef->getIndex() == index1node) {
				continue;
			}
			if((*it)->nodeRef->getIndex() == mstRoot->nodeRef->getIndex()) {
				continue;
			}
			
			double weight=(*it)->edge->getWeight();
			int source = (*it)->nodeRef->getIndex();
			
			if((*it)->parent==NULL) {
				continue;
			}
			
			int dest = (*it)->parent->nodeRef->getIndex();
			
			
		
			fb << (*it)->nodeRef->getIndex() <<" -- " << (*it)->parent->nodeRef->getIndex()<<
				"[label = \" "<< weight<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		
oneEdge1->getWeight()+oneEdge2->getWeight()
				+getPenalite(index1node, index1Edge1) + getPenalite(index1node, index1Edge2);
				
	
		double e1=oneEdge1->getWeight();
		fb << index1node <<" -- " << index1Edge1<<
				"[label = \" "<< e1<<"\" ]" <<	endl;
				
		double e2=oneEdge2->getWeight();
		fb << index1node <<" -- " << index1Edge2<<
				"[label = \" "<< e2<<"\" ]" <<	endl;

		fb << "}"<<endl;

        // Ferme le fichier :
        fb.close();
    }


}


void HeldKarp::printGrapheWithPenalties(string nomFichier) {
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"graph g{"<<endl;
    
      
    	
		for( graphEdges::const_iterator it = graphe->getEdges()->begin(); it!=graphe->getEdges()->end(); ++it) {
			//Edge ed = *it;
			double w=(*it)->getWeight() + getPenalite((*it)->getSource()->getIndex(), (*it)->getDest()->getIndex());
			if((*it)->isForced()) {
				w=std::numeric_limits<double>::infinity();
			}
			fb << (*it)->getSource()->getIndex() <<" -- " << (*it)->getDest()->getIndex()<<
				"[label = \" "<< w<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		fb << "}"<<endl;


        // Ferme le fichier :
        fb.close();
    }

}


std::vector<double> HeldKarp::getPenalites() {
	return penalite;
}
	
void HeldKarp::setPenalites(std::vector<double>& vecteur) {
	penalite.assign(vecteur.begin(), vecteur.end());
}



void HeldKarp::isTour(bool * isTour, bool * canContainTour) {
	
	*isTour = false;
	PrimMieux(canContainTour);
	
	if( (*canContainTour) == false)
		return;
	*isTour = true;
	double sum=0;
	for(std::vector<Pointeur*>::iterator it = pointeurs.begin(); it!= pointeurs.end() ; ++it) {
		int index = (*it)->nodeRef->getIndex();
		sum+=penalite.at(index);
		if((*it)->degre != 2) {
			*isTour = false;
		}

	}
	HKBound=oneTreeWeight-2*sum;
	return ;
}

