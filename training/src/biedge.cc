#include <binode.h> 


//Edge::Edge() {}

extern bool useReducedCostInAP;

BiEdge::BiEdge(BiNode * _leftNode,  BiNode * _rightNode, double _weight, Edge* _edge) :
	leftNode(_leftNode),
	rightNode(_rightNode),
	weight(_weight),
	edge(_edge),
	removed(true)
{
	if(leftNode->isInLeftSet() == false) {
		cout<<"Error. Try to add an edge with left node in right set."<<endl;
		exit(1);
	}
	if(rightNode->isInLeftSet() ) {
		cout<<"Error. Try to add an edge with right node in left set."<<endl;
		exit(1);
	}

	if(leftNode->getIndex() == rightNode->getIndex() ) {
		cout<<"Cannot add edge "<<print()<<endl;
		exit(1);
	}
}




/*
bool BiEdge::isTight() const{

//cout<<sourceNode->getPotential() + destNode->getPotential() <<endl;
//cout.flush();
//ATTENTION AUx ERREURS NUMERIQUES !!!!
	if( leftNode->getPotential() + rightNode->getPotential() <= weight + PRECISION 
	&&  leftNode->getPotential() + rightNode->getPotential() >= weight - PRECISION )
		return true;
		
	return false;
}
*/
double BiEdge::getReducedCost() const {
	double reducedCost = getWeight()-rightNode->getPotential()- leftNode->getPotential();
	if( reducedCost < -PRECISION) {
		cout<<"Found an edge with negative reduced cost. "<<reducedCost<<endl;
		cout<<print()<<endl;
		exit(1);
	}
	return reducedCost;
}



void BiEdge::setPositions(binodeEdges::iterator leftPos, binodeEdges::iterator rightPos, bigraphEdges::iterator graphPos) {
	leftPosition=leftPos;
	rightPosition=rightPos;
	graphPosition=graphPos;
}

binodeEdges::iterator BiEdge::getLeftPosition() const {
	return leftPosition;
}

binodeEdges::iterator BiEdge::getRightPosition() const {
	return rightPosition;
}

bigraphEdges::iterator BiEdge::getGraphPosition() const{
	return graphPosition;
}



BiNode * BiEdge::getRightNode() const{
	return rightNode;
}

BiNode * BiEdge::getLeftNode() const{
	return leftNode;
}

double BiEdge::getWeight() const {
	const double d=weight;
	if(useReducedCostInAP) {
		return edge->coutMarg;
	}
	else {
			return d;
	}
	return d;
}



Edge* BiEdge::getEdge() const {
	return edge;
}

std::string BiEdge::print() const{
/*
	cout<<"printEdge";
	cout.flush();
	*/
	stringstream uuh(stringstream::in | stringstream::out);
	uuh<<leftNode->print()<<" to "<<rightNode->print()<<"  "<<getWeight();
	
	/*
	cout<<"........done";
	cout.flush();
	
	*/
	return uuh.str();

}

bool BiEdge::isRemoved() {
	return removed;
}

void BiEdge::remove() {
	if(removed == true) {
		cout<<"Error. Removeing an edge twice from the graph."<<endl;
		exit(1);
	}
	removed=true;
}

void BiEdge::putBack() {
	if(removed == false) {
		cout<<"Error. Adding an edge twice to the graph."<<endl;
		exit(1);
	}
	removed=false;
}

