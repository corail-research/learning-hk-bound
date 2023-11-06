#include <node.h> 


//Edge::Edge() {}

Edge::Edge(Node* s, Node* d, double w) :
	source(s),
	dest(d),
	weight(w),
	forced(false),
	coutMarg(0.0),
	removed(true),
	biedge(NULL)
{


}


Edge::Edge(Edge& e) :
	source(e.source),
	dest(e.dest),
	weight(e.getWeight()),
	forced(e.isForced()),
	coutMarg(e.coutMarg)
{
};


void Edge::setPositions(nodeEdges::iterator sourcePos, nodeEdges::iterator destPos,graphEdges::iterator graphPos) {
	sourcePosition=sourcePos;
	destPosition=destPos;
	graphPosition=graphPos;
}

nodeEdges::iterator Edge::getSourcePosition() const {
	return sourcePosition;
}
nodeEdges::iterator Edge::getDestPosition() const{
	return destPosition;
}
graphEdges::iterator Edge::getGraphPosition() const{
	return graphPosition;
}

void Edge::remove() {
	if(removed == true) {
		cout<<"Error. Removeing an edge twice from the graph."<<endl;
		exit(1);
	}
	removed=true;
}

void Edge::putBack() {
	if(removed == false) {
		cout<<"Error. Adding an edge twice to the graph."<<endl;
		exit(1);
	}
	removed=false;
}

bool Edge::isRemoved() {
	return removed;
}
double Edge::getReplacementCost() {
	return replacement_cost;
}

void Edge::setReplacementCost(double cost) {
	replacement_cost=cost;
}

bool Edge::canBeRemoved() {
	if(source->getEdges()->size() <=2)
		return false;
	if(dest->getEdges()->size() <=2)
		return false;
		
	return true;
}

void Edge::force(bool * canContainTour) {
	if(forced) {
		cout<<"Error. Forcing an edge already forced."<<endl;
		exit(1);	
	}
	forced=true;

	source->addForcedVoisin();
	dest->addForcedVoisin();
	
	if(source->getNbForcedVoisins() >= 3 || dest->getNbForcedVoisins() >= 3 ) {
		*canContainTour = false;
		//cout<<"cannot force edge ";
		//print();
	}
	
}

void Edge::unforce() {
	if(!forced) {
		return;
	}
	forced=false;
	source->removeForcedVoisin();
	dest->removeForcedVoisin();
}

bool Edge::isForced() const {
	//cout <<forced<<endl;
	return forced;
}

Node * Edge::getSource() const{

	return source;
}

 Node * Edge::getDest() const{

	return dest;
}

double Edge::getWeight() {
	return weight;
}

double Edge::setWeight(double w) {
	weight = w;
}

BiEdge* Edge::getBiEdge() {
	return biedge;
}

void Edge::setBiEdge( BiEdge* _biedge) {
	
	biedge = _biedge;
}


bool  Edge::operator==(const Edge & b) const {
	return (*(this->source) == *(b.source)) && (*(this->dest) == *(b.dest));
}

bool  Edge::operator<(const Edge & b) const {
//on classe les edges par noeud source croissant,
//puis par noeud dest croissant
if(*(this->source) < *(b.source))
	return true;

if(*(b.source) < *(this->source) )
	return false;
	
//Sinon ils ont la mm source, on classe par dest
if(*(this->dest) < *(b.dest)  )
	return true;

//Si ils ont la meme dest, on renvoie false.
return false;

}

void Edge::print() const{
	//char* euh= new char[500];
	//sprintf(euh,"Edge from %d to %d  : %f \n", source.getIndex(), dest.getIndex(), weight);
	cout << "Edge from " << source->getIndex() << " to " << dest->getIndex() << " : "<< weight << endl;
	return ;

}
