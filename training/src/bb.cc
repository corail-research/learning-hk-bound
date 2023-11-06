#include <bb.h>
#include <limits>
#include <sstream>
#include <hungarianMethod.h>

extern bool heavyFiltering;
extern bool noSecondFiltering;
extern bool noFiltering;
extern bool forceCost;
extern bool forceCut;
extern bool forceDegree;
extern bool filterDegree;
extern bool printReducedGraph;
extern const char* reduceGraphFileName;

extern bool filterAP;
extern bool useAPBound;

extern bool verbose;

BB::BB(SymTSP* graph, double UB, int nF) :
graphe(graph),
hk(graph),
hg(graph),
upperBound(UB),
nbNode(0),
bb_tourFound(false),
tourCost(std::numeric_limits<double>::infinity()),
nbForce(nF),
nbForceCut(0),
nbForceCost(0),
percent_edges_filtered(-1.0)
{
	toForce_temp = new std::vector<Edge*>();
	toRemove_temp= new std::vector<Edge*>();
	
	toForce_temp->reserve(graph->getSize() );
	toRemove_temp->reserve(graph->getEdges()->size() - graph->getSize() );
	
	startTime =clock();

	
}
BB::~BB() {
	delete  toForce_temp;
	delete  toRemove_temp;
}
static bool filterMieux=false;

/*
void BB::test() {
		hk.index1node=1;

		hk.computeBound(upperBound);


		std::vector<Edge*>* toRemove = hk.filter(upperBound,filterMieux);
		graphe->removeEdges(toRemove);
		delete toRemove;
		
		
		std::vector<Edge*>* toForce=hk.forceMSTEdges(upperBound);
		graphe->forceEdges(toForce);
		delete toForce;
		
		
		pair<bool,std::vector<Edge*>* > pForceResult=hk.PrimForce();
		graphe->forceEdges(pForceResult.second);
		cout<<pForceResult.second->size()<<" arc forces apres Prim force"<<endl;
		delete pForceResult.second;


		hk.PrimMieux();
		
		

		
		toRemove = hk.filter(upperBound,filterMieux);
		graphe->removeEdges(toRemove);
		delete toRemove;

		toForce=hk.forceMSTEdges(upperBound);
		graphe->forceEdges(toForce);
		delete toForce;
		

		
		pForceResult=hk.PrimForce();
		graphe->forceEdges(pForceResult.second);
		cout<<pForceResult.second->size()<<" arc forces apres Prim force"<<endl;
		delete pForceResult.second;

		
		
		
		

		toRemove = hk.filter(upperBound,filterMieux);
		graphe->removeEdges(toRemove);
		delete toRemove;
		
		toForce=hk.forceMSTEdges(upperBound);
		graphe->forceEdges(toForce);
		delete toForce;
		
		pForceResult=hk.PrimForce();
		graphe->forceEdges(pForceResult.second);
		cout<<pForceResult.second->size()<<" arc forces apres Prim force"<<endl;
		delete pForceResult.second;
		
		
		
		
		toRemove = hk.filter(upperBound,filterMieux);
		graphe->removeEdges(toRemove);
		delete toRemove;
		
		toForce=hk.forceMSTEdges(upperBound);
		graphe->forceEdges(toForce);
		delete toForce;
		
		pForceResult=hk.PrimForce();
		graphe->forceEdges(pForceResult.second);
		cout<<pForceResult.second->size()<<" arc forces apres Prim force"<<endl;
		delete pForceResult.second;



	//	hk.printDot("euh.dot");

	cout <<"Remaining edges : "<<graphe->getEdges()->size()<<endl;
	cout<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;
}	

*/

void BB::compute(int size_bf, float bound_factors[]) {
	
	int initial_nb_edges = graphe->getEdges()->size();

	double maxBound=0.0;
	int oneNodeMaxBound=0;
	//Premier filtrage
	hk.index1node=0;
	
	std::vector<Edge*>* toRemove = new std::vector<Edge*>();
	std::vector<Edge*>* toForce = new std::vector<Edge*>();

	float HKbound_fnf; 
	float factors[size_bf-1];
	bound_factors[0] = -1;

	bool canStopBranching = false;
	
	//try hk with all nodes as one-node, and keep the node giving the best bound
	for(int i=0; i<1; i++) { //graphe->getSize()
		if(verbose) cout<<"One Node :"<<i<<endl;
		hk.index1node=i;
		
		//filter_and_force(toRemove, toForce,  &canStopBranching);
		filter_and_force(toRemove, toForce,  &canStopBranching, &HKbound_fnf, factors);
		
		if(HKbound_fnf > bound_factors[0]){
			bound_factors[0] = HKbound_fnf;
			for(int j=0; j<size_bf; j++){
				bound_factors[j+1] = factors[j];
			}
		}
		////if HKbound_fnf > bound_factors[0]: bound_factors = [HKbound_fnf,factors]
		
		toRemove->clear();
		toForce->clear();
	
		
		if(canStopBranching) {
			if(verbose) cout<<"can stop branching"<<endl;
			break;
		}
		
		/*if(hk.getBound() > maxBound) {
			maxBound=hk.getBound();
			oneNodeMaxBound=i;
		}*/
		
	}



	//RETURN BOUND AND PENALTIES

	
	/*int current_nb_edges = graphe->getEdges()->size();
	percent_edges_filtered = 100.0*(initial_nb_edges-current_nb_edges)/(1.0 * initial_nb_edges);
	
	if(canStopBranching == false) {

		if(verbose) cout <<endl<<endl<<"Best borne :"<<maxBound<<"   avec one node "<<oneNodeMaxBound<< endl;
		hk.index1node=oneNodeMaxBound;
		
		filter_and_force(toForce, toRemove,  &canStopBranching);
		
		delete toForce;
		delete toRemove;

	
	
		dfs_Remove(0);
	}*/
	
	
	/*
	if(printReducedGraph) {
		graphe->printInTSPLIBFormat(reduceGraphFileName);
		cout<<"Reduced graph printed in "<<reduceGraphFileName<<endl;
		exit(0);
	}
		*/
	

	

	
	printEndInfos(true);
	
	/*
	if(bb_tourFound) {
			cout<<endl<<endl<<" Tour found :"<<tourCost<<endl;
			cout<<"UpperBound : "<<upperBound<<endl;
			
			cerr<<" Tour found :"<<tourCost<<endl;
			cerr<<"UpperBound : "<<upperBound<<endl;
	}
	else {
		cout<<endl<<endl<<"No tour found :("<<endl;
		
		cerr<<"No tour found :("<<endl;
	}
	*/
	cout<<endl<<"Nb BB nodes : "<<nbNode<<endl;
	cerr<<"Nb BB nodes : "<<nbNode<<endl<<endl<<endl;
	
	
	//cerr<<"Nb force cut "<<nbForceCut<<endl
	//<<"Nb force cost "<<nbForceCost<<endl;
}





bool coutMargComparator (Edge* u, Edge* v) {
	/*
	if(u->getSource()->getEdges()->size() <3 || u->getDest()->getEdges()->size() <3 ){
		if(v->getSource()->getEdges()->size() <3 || v->getDest()->getEdges()->size() <3 ) {
			return u->coutMarg > v->coutMarg;
		}
		return true;
	}
	if(v->getSource()->getEdges()->size() <3 || v->getDest()->getEdges()->size() <3 ){
		return false;
	}
	
	*/
	
	
	return u->coutMarg > v->coutMarg;
	
	/*
	int a =u->getSource()->getEdges()->size();
	if(u->getDest()->getEdges()->size() < a)
		a=u->getDest()->getEdges()->size();
		
	int b =v->getSource()->getEdges()->size();
	if(v->getDest()->getEdges()->size() < b)
		b=v->getDest()->getEdges()->size();
		
	return a< b;
	*/
	
}



bool replacementCostComparator (Edge* u, Edge* v) {

	
	return u->getReplacementCost() > v->getReplacementCost();
	

	
}

pair<bool,Edge*> BB::getEdgeWithMaxReplacementCost() {

	//pair<bool,bool> hkResult=hk.computeBound(upperBound);
	/*
	bool tourFound=hk.isTour();
	
	if(tourFound) {
		cout<<"Found a tour of size "<<hk.getBound()<<endl;
		
		if(hk.getBound() <tourCost) {
			tourCost=hk.getBound();
			if(tourCost <upperBound) {
				upperBound=tourCost;
			}
			//hk.printDot("tour.dot");
		}
		//return;
		//exit(0);
	}
	*/
	
	//if no filtering was done, we need to filter in order to compute replacement costs
	//if APbound was used, replacement costs has changed ??
	/*
	if(noFiltering || useAPBound) {
		bool canStopBranching = false;
		hk.PrimMieux(&canStopBranching);
		if( canStopBranching == false) {
			cout<<"Looking for edge with max replacement cost whereas can stop branching."<<endl;
			cout<<"Should have been detected sooner."<<endl;
			exit(1);
		}
		toRemove_temp->clear();
		hk.filter(upperBound,toRemove_temp);
		toRemove_temp->clear();
	}
	* */
	bool edgeFound=false;
	Edge* maxEdge;
	double maxCost=-std::numeric_limits<double>::infinity();


	
	for(std::vector<Pointeur*>::iterator it = hk.pointeurs.begin(); it!= hk.pointeurs.end(); ++it) {
		if((*it)->rank == 0 ){ //on ne regarde pas la racine de l'arbre
			continue;
		}
	
		if((*it)->nodeRef->getIndex() == hk.index1node ){ //on ne regarde pas le oneNode
			continue;
		}
		if((*it)->edge== NULL ) {
			cout<<"erreur, node sans parent dans le mst, c'est mauvais"<<endl;
			cout<<"node "<<(*it)->nodeRef->getIndex()<<endl;
			cout<<"rang "<<(*it)->rank<<endl;
			
			hk.printDot("uh.dot");
			exit(1);
		}
		if((*it)->edge->isForced() ) {//on ne regarde pas les arcs deja forcé
			continue;
		}
		
		double replacementCost=(*it)->minEdgeWeight - (*it)->realWeight;
		
		if(replacementCost > maxCost) {
			maxCost=replacementCost;
			maxEdge=(*it)->edge;
			edgeFound=true;
		}

	
	}
	
	//faut ajouter oneEdge1 et oneEdge2 !
	//Si il existe une troisieme arete connectée au one node 
	if(hk.oneEdge3Weight  <  std::numeric_limits<double>::infinity() ) {
	
		if(!hk.oneEdge1->isForced() ) {
		edgeFound=true;
		double replacementCost=hk.oneEdge3Weight - hk.oneEdge1Weight;
		if(replacementCost > maxCost) {
			maxCost=replacementCost;
			maxEdge=hk.oneEdge1;
		}

		}
		if(!hk.oneEdge2->isForced() ) {
		edgeFound=true;
		double replacementCost=hk.oneEdge3Weight - hk.oneEdge2Weight;
		if(replacementCost > maxCost) {
			maxCost=replacementCost;
			maxEdge=hk.oneEdge2;
		}
		}
	}
	

	return pair<bool,Edge*>(edgeFound,maxEdge);

}



void BB::getEdgesToBranchOn(std::vector<Edge*>* edgesToBranchOn) {

	//pair<bool,bool> hkResult=hk.computeBound(upperBound);
	/*
	bool tourFound=hk.isTour();
	
	if(tourFound) {
		cout<<"Found a tour of size "<<hk.getBound()<<endl;
		
		if(hk.getBound() <tourCost) {
			tourCost=hk.getBound();
			if(tourCost <upperBound) {
				upperBound=tourCost;
			}
			//hk.printDot("tour.dot");
		}
		//return;
		//exit(0);
	}
	*/
	
	
	
	//if no filtering was done, we need to filter in order to compute replacement costs
	//if APbound was used, replacement costs has changed ??
	/*
	bool canStopBranching = false;
	hk.PrimMieux(&canStopBranching);
	if( canStopBranching == false) {
		cout<<"Looking for edge with max replacement cost whereas can stop branching."<<endl;
		cout<<"Should have been detected sooner."<<endl;
		exit(1);
	}
	toRemove_temp->clear();
	hk.filter(upperBound,toRemove_temp);
	toRemove_temp->clear();
	
*/

	
	for(std::vector<Pointeur*>::iterator it = hk.pointeurs.begin(); it!= hk.pointeurs.end(); ++it) {
		if((*it)->rank == 0 ){ //on ne regarde pas la racine de l'arbre
			continue;
		}
	
		if((*it)->nodeRef->getIndex() == hk.index1node ){ //on ne regarde pas le oneNode
			continue;
		}
		if((*it)->edge== NULL ) {
			cout<<"erreur, node sans parent dans le mst, c'est mauvais"<<endl;
			cout<<"node "<<(*it)->nodeRef->getIndex()<<endl;
			cout<<"rang "<<(*it)->rank<<endl;
			
			//hk.printDot("uh.dot");
			exit(1);
		}
		if((*it)->edge->isForced() ) {//on ne regarde pas les arcs deja forcé
			continue;
		}
		
	
		(*it)->edge->setReplacementCost((*it)->minEdgeWeight - (*it)->realWeight);
		
		edgesToBranchOn->push_back((*it)->edge);
		
		if((*it)->edge->isRemoved() ) {
			cout<<"There is a removed edge in the onetree"<<endl;
			(*it)->edge->print();
			exit(1);
		}
	
	}
	
	//faut ajouter oneEdge1 et oneEdge2 !
	//Si il existe une troisieme arete connectée au one node 
	if(hk.oneEdge3Weight  <  std::numeric_limits<double>::infinity() ) {
	
		if(!hk.oneEdge1->isForced() ) {
			hk.oneEdge1->setReplacementCost(hk.oneEdge3Weight - hk.oneEdge1Weight);
			edgesToBranchOn->push_back(hk.oneEdge1);
			
		}
		if(!hk.oneEdge2->isForced() ) {
			hk.oneEdge2->setReplacementCost(hk.oneEdge3Weight - hk.oneEdge2Weight);
			edgesToBranchOn->push_back(hk.oneEdge2);
		}
	}
	

	return;

}





void BB::dfs_Remove(int profondeur) {
	
	
	nbNode++;
	
	
	
	if(verbose) 	cout<<endl<<"Node "<<nbNode<<endl;
	int thisNode=nbNode;
	if(verbose) cout<<"Profondeur dfs : "<<profondeur<<endl;
	if(verbose) cout<<"Nb arcs restants : "<<graphe->getEdges()->size()<<endl;
	if(verbose) cout<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;



	std::vector<Edge*>* toRemove =new std::vector<Edge*>();
	std::vector<Edge*>* toForce =new std::vector<Edge*>();
	
	
	bool canStopBranching = false;
	
	//hk.clearPenalties();
	
	//filter_and_force(toRemove, toForce, &canStopBranching);
	

	std::vector<double> penalites =hk.getPenalites();
	
	if(canStopBranching) {
		graphe->addEdges(toRemove);
		hg.addEdges(toRemove);
		delete toRemove;
		
		graphe->unforceEdges(toForce);
		delete toForce;	
		
		//cout<<endl<<"Node "<<nbNode<<"  can stop branching"<<endl;
		return;
	}
	
	
	
		
	//continuer la dfs
	std::vector<Edge*>* edgesToBranchOn = new std::vector<Edge*>();
	//used computed reduced costs to sort edges to branch on
	//getEdgesToBranchOn(edgesToBranchOn);
	//sort edges based on replacement costs
	//std::sort(edgesToBranchOn->begin(), edgesToBranchOn->end(), replacementCostComparator );
	

	int nbVu=0;
	
	pair<bool,Edge*> nextEdge=getEdgeWithMaxReplacementCost();
	
	if(nextEdge.first) {
	

		Edge* e=nextEdge.second;
		
		nbVu++;
		if(verbose) cout<<"Node "<<thisNode<<"  "<<nbVu<<endl;
		

		bool canContainTour = true;
		graphe->removeEdge(e, &canContainTour);
		hg.removeEdge(e);

		//if removing this edge don't create a node of degree 1 or 0,
		//remove this edge and branch
		if(canContainTour) {
			if(verbose) cout<<"BB-Remove edge :";
			if(verbose) e->print();
			
			dfs_Remove(profondeur+1);
			
		//	hk.setPenalites(penalites);
			//hk.clearPenalties();
		}
		else {
			if(verbose) cout<<"BB-Cannot remove edge ";
			if(verbose) 	e->print();
			if(verbose) cout<<"  "<<e->getSource()->getEdges()->size()<<"  "<<e->getDest()->getEdges()->size()<<endl;
		}
		
		graphe->addEdge(e);
		hg.addEdge(e);
		canContainTour = true;
		e->force(&canContainTour);
		toForce->push_back(e);
		
		//if forcing this edge  create a node with 3 forced edge, stop branching
		if(canContainTour == false) {
			if(verbose) cout<<"BB-Cannot force edge ";
			if(verbose) e->print();
			//break;
		}
		else {
			dfs_Remove(profondeur);
		}
		 
	}
	if(verbose) cout<<"Node "<<thisNode<<"end"<<nbVu<<endl;
	graphe->addEdges(toRemove);
	hg.addEdges(toRemove);
	delete toRemove;

	graphe->unforceEdges(toForce);
	delete toForce;
	
	delete edgesToBranchOn;


	return;

}

void BB::filter_and_force(std::vector<Edge*>* toRemove, std::vector<Edge*>* toForce, bool *canStopBranching, float *HKbound_fnf, float factors[]) {
	clock_t endTime=clock();

	int mn =(int) ( (endTime-startTime)/(60.0*(double)CLOCKS_PER_SEC) );
	
	if(mn > 30) {
		cerr<<"Stopping after "<<mn<<" mn."<<endl;
		printEndInfos(false);
		
		//cerr<<"Nb BB nodes : "<<nbNode<<endl;
		//exit(1);
	}
	
	*canStopBranching = false;
	bool tourFound = false;
	bool canContainTour = true;
	
	//hk bound
	//hk.computeBound(upperBound, &tourFound, canStopBranching);
	hk.computeBound(upperBound, &tourFound, canStopBranching, HKbound_fnf, factors);

		
	/*

	//if a tour was found,
	//update upper bound and stop branching
	if(tourFound) {
		bb_tourFound= true;
		if(verbose) cout<<"Found a tour of size "<<hk.getBound()<<endl;
		if(hk.getBound() <tourCost) {
			tourCost=hk.getBound();
			if(tourCost <upperBound) {
				upperBound=tourCost;
			}
		}
		*canStopBranching = true;
		return ;
	}
	
	if( (*canStopBranching) == true) {
		return;
	}
	
	
	
	

	
	
	//compute reduce costs and replacement costs
	toRemove_temp->clear();
	hk.filter_reducedCosts(upperBound,toRemove_temp);
	
	//filter edges based on reduced costs if asked
	if(noFiltering == false) {
		graphe->removeEdges(toRemove_temp, &canContainTour);
		hg.removeEdges(toRemove_temp);
		toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
		if( canContainTour == false) {
			*canStopBranching = true;
			return;
		}
	}
	toRemove_temp->clear();
	
	//if nbForce==-1, we repeat until fixed point is reached

	bool hasChanged=true;
	int nbTours=0;
	while(hasChanged) {
		nbTours++;
		hasChanged=false;
	
		//force edges based on replacement costs
		if(forceCost) {
			if(verbose) cout<<"Force with cost method.  ";
			toForce_temp->clear();
			hk.forceMSTEdges(upperBound, toForce_temp);
			graphe->forceEdges(toForce_temp, &canContainTour);
			toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
			if( !toForce_temp->empty() ) {
				hasChanged=true;
			}
			toForce_temp->clear();
			
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
		}

		
		
		//compute AP bound if asked
		if( graphe->isFromATSP && useAPBound) {
			hg.computeBound(&canContainTour);
			
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
			double addBound=hg.getBound();
			if(verbose) cout<<"AP bound : "<<addBound<<"   Total bound :  "<<hk.getBound()+addBound<<endl;
			
			if(hk.getBound()+addBound > upperBound + PRECISION ) {
				*canStopBranching = true;
				return;
			}
		}
		
		
		
		//filter edges based on AP bound
		if( graphe->isFromATSP && filterAP) {
			toRemove_temp->clear();
			hg.filter(upperBound, hk.getBound(), toRemove_temp);
			graphe->removeEdges(toRemove_temp, &canContainTour);
			hg.removeEdges(toRemove_temp);
			toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
			if( !toRemove_temp->empty() ) {
				hasChanged=true;
			}
			toRemove_temp->clear();
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
			
		}
		
	
		//force edges based on cuts' size
		if(forceCut) {
			if(verbose) cout<<"Force with cut method.  ";
			toForce_temp->clear();
			hk.PrimForce(&canContainTour, toForce_temp);
			graphe->forceEdges(toForce_temp, &canContainTour);
			toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
			if( !toForce_temp->empty() ) {
				hasChanged=true;
			}
			toForce_temp->clear();
			
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}

		}



		//force edges based on nodes' degree
		if(forceDegree) {
			if(verbose) cout<<"Force with degree method.  ";
			toForce_temp->clear();
			hk.forceEdgesWithDegree2(toForce_temp);
			graphe->forceEdges(toForce_temp, &canContainTour);
			toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
			if( !toForce_temp->empty() ) {
				hasChanged=true;
			}
			toForce_temp->clear();
			
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
		}
	
		//filter edges based on forced degree
		if(filterDegree) {
			toRemove_temp->clear();
			hk.filter_forcedDegrees(toRemove_temp);
			graphe->removeEdges(toRemove_temp, &canContainTour);
			hg.removeEdges(toRemove_temp);
			toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
			if( !toRemove_temp->empty() ) {
				hasChanged=true;
			}
			toRemove_temp->clear();
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
			
		}

		//forcing edges with cut and degree method can force
		//edges which were not in the one-tree
		//filtering edges with AP and filterDegree methods
		//can filter edges which were inside the one-tree
		//thus the one-tree is recomputed 
		hk.isTour(&tourFound, &canContainTour);
		
		if(tourFound) {
			if(verbose) cout<<"Found a tour of size "<<hk.getBound()<<endl;
	
			bb_tourFound=true;
			if(hk.getBound() <tourCost) {
				tourCost=hk.getBound();
				if(tourCost <upperBound) {
					upperBound=tourCost;
				}
			}
			*canStopBranching = true;
			return;
		}
		if( canContainTour == false) {
			*canStopBranching = true;
			return;
		}
			
		
		
		//filter again	
		
		//recompute reduced costs and replacement costs
		toRemove_temp->clear();
		hk.filter_reducedCosts(upperBound,toRemove_temp);
		//remove edges based on reduced costs if asked
		if( noFiltering == false ) {
			graphe->removeEdges(toRemove_temp, &canContainTour);
			hg.removeEdges(toRemove_temp);
			toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
			if( !toRemove_temp->empty() ) {
				hasChanged=true;
			}
			toRemove_temp->clear();
			if( canContainTour == false) {
				*canStopBranching = true;
				return;
			}
		}
		else {
			toRemove_temp->clear();
		}
		
		if(nbForce > 0  && nbTours >= nbForce )	{
			break;
		}
	}

	
	if(verbose) cout<<"Remaining edges : "<<graphe->getEdges()->size()<<endl;
	if(verbose) cout<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;
	
	//cerr<<"First bound : "<<hk.getBound()<<endl;
	//cerr<<"UB : "<<upperBound<<endl;
	//cerr<<"Remaining edges : "<<graphe->getEdges()->size()<<endl;
	//cerr<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;
	
	
	*/
	
}




void BB::compute_ap() {
	int initial_nb_edges = graphe->getEdges()->size();
	if( graphe->isFromATSP == false ) {
		cout<<"AP bound works only with atsp."<<endl;
		exit(1);
	}

	double maxBound=0.0;
	int oneNodeMaxBound=0;
	//Premier filtrage
	hk.index1node=0;
	
	std::vector<Edge*>* toRemove = new std::vector<Edge*>();
	std::vector<Edge*>* toForce = new std::vector<Edge*>();

	bool canStopBranching = false;
	

	hk.index1node=0;
	
	filter_and_force_ap(toRemove, toForce,  &canStopBranching);
	
	toRemove->clear();
	toForce->clear();



	filter_and_force_ap(toForce, toRemove,  &canStopBranching);
	
	delete toForce;
	delete toRemove;

	int current_nb_edges = graphe->getEdges()->size();
	percent_edges_filtered = 100.0* current_nb_edges/(1.0 * initial_nb_edges);
	

	
	if(canStopBranching == false) {
		dfs_Remove_ap(0);
	}
	
	
	/*
	if(printReducedGraph) {
		graphe->printInTSPLIBFormat(reduceGraphFileName);
		cout<<"Reduced graph printed in "<<reduceGraphFileName<<endl;
		exit(0);
	}
		*/
	

	

	
	printEndInfos(true);
	/*
	if(bb_tourFound) {
			cout<<endl<<endl<<" Tour found :"<<tourCost<<endl;
			cout<<"UpperBound : "<<upperBound<<endl;
			
			cerr<<" Tour found :"<<tourCost<<endl;
			cerr<<"UpperBound : "<<upperBound<<endl;
	}
	else {
		cout<<endl<<endl<<"No tour found :("<<endl;
		
		cerr<<"No tour found :("<<endl;
	}
	
	cout<<endl<<"Nb BB nodes : "<<nbNode<<endl;
	cerr<<"Nb BB nodes : "<<nbNode<<endl<<endl<<endl;
	*/
	
	//cerr<<"Nb force cut "<<nbForceCut<<endl
	//<<"Nb force cost "<<nbForceCost<<endl;
}



void BB::filter_and_force_ap(std::vector<Edge*>* toRemove, std::vector<Edge*>* toForce, bool *canStopBranching) {
	
	clock_t endTime=clock();

	int mn =(int) ( (endTime-startTime)/(60.0*CLOCKS_PER_SEC) );
	
	if(mn >= 30) {
		cerr<<"Stopping after "<<mn<<" mn."<<endl;
		printEndInfos(false);
		//cerr<<"Nb BB nodes : "<<nbNode<<endl;
		exit(1);
	}
	*canStopBranching = false;
	bool tourFound = false;
	bool canContainTour = true;
	
	bool hasChanged = false;
	
	//force edges based on cuts' size
	if(forceCut) {
		if(verbose) cout<<"Force with cut method.  ";
		toForce_temp->clear();
		hk.PrimForce(&canContainTour, toForce_temp);
		graphe->forceEdges(toForce_temp, &canContainTour);
		toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
		toForce_temp->clear();
		
		if( canContainTour == false) {
			if(verbose) cout<<"force cut->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}

	}
	
	//force edges based on nodes' degree
	if(forceDegree) {
		if(verbose) cout<<"Force with degree method.  ";
		toForce_temp->clear();
		hk.forceEdgesWithDegree2(toForce_temp);
		graphe->forceEdges(toForce_temp, &canContainTour);
		toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
		toForce_temp->clear();
		
		if( canContainTour == false) {
			if(verbose) cout<<"force degree->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
	}

	//filter edges based on forced degree
	if(filterDegree) {
		if(verbose) cout<<"Filter with degree method.  ";
		toRemove_temp->clear();
		hk.filter_forcedDegrees(toRemove_temp);
		//cout<<toRemove_temp->size()<<" "<<endl;
		graphe->removeEdges(toRemove_temp, &canContainTour);
		hg.removeEdges(toRemove_temp);
		toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
		toRemove_temp->clear();
		if( canContainTour == false) {
			if(verbose) cout<<"filter degree->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
		
	}
	
	//compute AP bound 
	hg.computeBound(&canContainTour);
	
	if( canContainTour == false) {
		if(verbose) cout<<"AP->can stop branching: "<<endl;
		*canStopBranching = true;
		return;
	}
	double addBound=hg.getBound();
	if(verbose) cout<<"AP bound : "<<addBound<<endl;
	if(hg.isTour() == true) {
			bb_tourFound = true;
			if(addBound < tourCost) {
				tourCost = addBound;
				if(tourCost <upperBound) {
					upperBound=tourCost;
				}
			}
	}
	
	if(addBound > upperBound + PRECISION ) {
		if(verbose) cout<<"AP bound > UB: "<<endl;
		*canStopBranching = true;
		return;
	}
	
	
	//filter edges based on AP bound
	if(filterAP) {
		if(verbose) cout<<"Filter AP.  ";
		toRemove_temp->clear();
		hg.filter(upperBound, 0.0, toRemove_temp);
		graphe->removeEdges(toRemove_temp, &canContainTour);
		hg.removeEdges(toRemove_temp);
		toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
		toRemove_temp->clear();
		if( canContainTour == false) {
			if(verbose) cout<<"filter AP->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
		
	}
	
	
	//set potentials in hk and compute a one-tree
	hg.setEdgesCostsInHK();
	hk.isTour(&tourFound, &canContainTour);
	
	if(tourFound) {
		if(verbose) cout<<"Found a tour of size "<<hk.getBound()<<endl;

		bb_tourFound=true;
		if(hk.getBound() + addBound <tourCost) {
			tourCost=hk.getBound()  + addBound;
			if(tourCost <upperBound) {
				upperBound=tourCost;
			}
		}
		*canStopBranching = true;
		return;
	}
	if( canContainTour == false) {
		if(verbose) cout<<"hk->can stop branching: "<<endl;
		*canStopBranching = true;
		return;
	}
	if(verbose) cout<<"HK bound "<<hk.getBound()<<" Total bound "<<hk.getBound() +addBound<<endl;
	if(hk.getBound() +addBound> upperBound + PRECISION) {
		if(verbose) cout<<"hk-> LB > UB: "<<endl;
		*canStopBranching = true;
		return;
	}
	
	//compute reduced costs and replacement costs
	if( noFiltering == false || forceCost) {
		toRemove_temp->clear();
		hk.filter_reducedCosts(upperBound-addBound,toRemove_temp);
	}
	//remove edges based on reduced costs if asked
	if( noFiltering == false ) {
		if(verbose) cout<<"Filter cost.  ";
		graphe->removeEdges(toRemove_temp, &canContainTour);
		hg.removeEdges(toRemove_temp);
		toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
		if(toRemove_temp->size() > 0 ){
			hasChanged = true;
		}
		toRemove_temp->clear();
		if( canContainTour == false) {
			if(verbose) cout<<"filter cost->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
	}
	else {
		toRemove_temp->clear();
	}
	
	
	//force edges based on replacement costs
	if(forceCost) {
		if(verbose) cout<<"Force with cost method.  ";
		toForce_temp->clear();
		hk.forceMSTEdges(upperBound-addBound, toForce_temp);
		graphe->forceEdges(toForce_temp, &canContainTour);
		toForce->insert(toForce->end(), toForce_temp->begin(), toForce_temp->end() );
		if(toForce_temp->size() > 0 ){
			hasChanged = true;
		}
		toForce_temp->clear();
		
		if( canContainTour == false) {
			if(verbose) cout<<"force cost->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
	}
	

	
	


	

	
	if(hasChanged) {
		//recompute AP bound 
		hg.computeBound(&canContainTour);
		
		if( canContainTour == false) {
			if(verbose) cout<<"AP->can stop branching: "<<endl;
			*canStopBranching = true;
			return;
		}
		addBound=hg.getBound();
		if(verbose) cout<<"AP bound : "<<addBound<<endl;;
		
		if(addBound > upperBound + PRECISION ) {
			*canStopBranching = true;
			return;
		}
		
		
		//refilter edges based on AP bound
		if(filterAP) {
			if(verbose) cout<<"Filter AP.  ";
			toRemove_temp->clear();
			hg.filter(upperBound, 0.0, toRemove_temp);
			graphe->removeEdges(toRemove_temp, &canContainTour);
			hg.removeEdges(toRemove_temp);
			toRemove->insert(toRemove->end(), toRemove_temp->begin(), toRemove_temp->end() );
			toRemove_temp->clear();
			if( canContainTour == false) {
				if(verbose) cout<<"filter AP->can stop branching: "<<endl;
				*canStopBranching = true;
				return;
			}
			
		}
	
	}

	
	

	
	if(verbose) cout<<"Remaining edges : "<<graphe->getEdges()->size()<<endl;
	if(verbose) cout<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;
	
	//cerr<<"First bound : "<<hk.getBound()<<endl;
	//cerr<<"UB : "<<upperBound<<endl;
	//cerr<<"Remaining edges : "<<graphe->getEdges()->size()<<endl;
	//cerr<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;
	
	
	
	
}



void BB::dfs_Remove_ap(int profondeur) {
	
	
	nbNode++;
	
	if(verbose) cout<<endl<<"Node "<<nbNode<<endl;
	int thisNode=nbNode;
	if(verbose) cout<<"Profondeur dfs : "<<profondeur<<endl;
	if(verbose) cout<<"Nb arcs restants : "<<graphe->getEdges()->size()<<endl;
	if(verbose) cout<<"Nb arcs forces : "<<graphe->getNbArcForce()<<endl;



	std::vector<Edge*>* toRemove =new std::vector<Edge*>();
	std::vector<Edge*>* toForce =new std::vector<Edge*>();
	
	
	bool canStopBranching = false;
	
	//hk.clearPenalties();
	
	filter_and_force_ap(toRemove, toForce, &canStopBranching);
	

	if(canStopBranching) {
		if(verbose) 	cout<<"Node "<<thisNode<<" end "<<endl;
		graphe->addEdges(toRemove);
		hg.addEdges(toRemove);
		delete toRemove;
		
		graphe->unforceEdges(toForce);
		delete toForce;	
		
		//cout<<endl<<"Node "<<nbNode<<"  can stop branching"<<endl;
		return;
	}
	
	
	
		
	//continuer la dfs
	std::vector<Edge*>* edgesToBranchOn = new std::vector<Edge*>();
	//used computed reduced costs to sort edges to branch on
	hg.getEdgesToBranchOn(edgesToBranchOn);
	//sort edges based on replacement costs
	std::sort(edgesToBranchOn->begin(), edgesToBranchOn->end(), replacementCostComparator );
	

	int nbVu=0;

	
	
	Edge* e=(*(edgesToBranchOn->begin()));
	
	nbVu++;
	if(verbose) cout<<"Node "<<thisNode<<"  "<<nbVu<<endl;
	

	bool canContainTour = true;
	graphe->removeEdge(e, &canContainTour);
	hg.removeEdge(e);

	//if removing this edge don't create a node of degree 1 or 0,
	//remove this edge and branch
	if(canContainTour) {
		if(verbose) cout<<"BB-Remove edge :";
		if(verbose) e->print();
		
		dfs_Remove_ap(profondeur+1);
		
	}
	else {
		if(verbose) cout<<"BB-Cannot remove edge ";
		if(verbose) e->print();
		if(verbose) cout<<"  "<<e->getSource()->getEdges()->size()<<"  "<<e->getDest()->getEdges()->size()<<endl;
	}
	
	graphe->addEdge(e);
	hg.addEdge(e);
	canContainTour = true;
	e->force(&canContainTour);
	toForce->push_back(e);
	
	//if forcing this edge  create a node with 3 forced edge, stop branching
	if(canContainTour == false) {
		if(verbose) cout<<"BB-Cannot force edge ";
		if(verbose) e->print();
		
	}
	else {
		dfs_Remove_ap(profondeur+1);
	}
		 
	
	if(verbose) cout<<"Node "<<thisNode<<" end "<<nbVu<<endl;
	graphe->addEdges(toRemove);
	hg.addEdges(toRemove);
	delete toRemove;

	graphe->unforceEdges(toForce);
	delete toForce;
	
	delete edgesToBranchOn;


	return;

}

void BB::printEndInfos(bool isOptimal) {
	if(bb_tourFound) {
			cout<<endl<<endl<<" Tour found : "<<tourCost<<endl;
			cout<<"UpperBound : "<<upperBound<<endl;
			
			cerr<<" Tour found : "<<tourCost<<endl;
			cerr<<"UpperBound : "<<upperBound<<endl;
	}
	else {
		cout<<endl<<endl<<" Tour found : NO"<<endl;
		
		cerr<<" Tour found : NO"<<endl;
	}
	
	cout<<"Optimality proved : "<<isOptimal<<endl;
	cerr<<"Optimality proved : "<<isOptimal<<endl;
	
	cout<<endl<<"Nb BB nodes : "<<nbNode<<endl;
	cerr<<"Nb BB nodes : "<<nbNode<<endl<<endl<<endl;
	
	cerr<<"percent edges filtered before branching : "<<percent_edges_filtered<<endl;
	
	//exit(0);
}
