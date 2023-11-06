#include <binode.h> 



BiNode::BiNode(int _index,bool _isInLeftSet,double _potential, Node* _masterNode) :
	index(_index),
	masterNode(_masterNode),
	set(_isInLeftSet),
	potential(_potential),
	slack(std::numeric_limits<double>::infinity()),
	parentEdge(NULL),
	matchingEdge(NULL),
	inTree(false),
	position_listNotInTree(NULL)
{
	edges = new binodeEdges();

};

BiNode::~BiNode() {
		delete edges;
}

binodeEdges* BiNode::getEdges() {
	return  edges;
}



binodeEdges::iterator BiNode::addEdge( BiEdge * e) {
	//Check if edge e is connected to this node 
	if(e->getLeftNode()->getIndex()!= index && e->getRightNode()->getIndex()!= index) {
		cout<<"Cannot add edge "<<e->print()<<" to node "<<index<<endl;
		exit(1);
	}

	edges->push_back(e);	
	return --edges->end();
}



void BiNode::removeEdge(binodeEdges::iterator position){
	edges->erase(position );

}


double BiNode::getPotential() const{
	return potential;
}

void BiNode::setPotential(double _potential) {
	potential=_potential;
}


int BiNode::getIndex() const {
	int i=index;
	return i;
}

bool BiNode::isInLeftSet() const {
	return set;
}

Node * BiNode::getMasterNode() {
	Node * result = masterNode;
	return result;
}


BiEdge * BiNode::getParentEdge() {
		return parentEdge;
}
/*
BiNode * BiNode::getNodeMatchedTo() {
		return matchedTo;
}
*/
void BiNode::setParentEdge(BiEdge * edge) {
	if(set) {
		if(edge->getLeftNode()->getIndex() != index) {
			cout<<"Error. Try to set parent of node "<<index<<" using edge "<<edge->print()<<endl;
			exit(1);
		}
	}
	else {
		if(edge->getRightNode()->getIndex() != index) {
			cout<<"Error. Try to set parent of node "<<index<<" using edge "<<edge->print()<<endl;
			exit(1);
		}
	}
	parentEdge = edge;
}
/*
void BiNode::matchTo(BiNode * node) {
	if(node->isInLeftSet() == isInLeftSet() ) {
		cout<<"Try to match two nodes in the same set"<<endl;
		exit(1);
	}	
	matchedTo = node;
}
*/
void BiNode::clear() {
	//cout<<"clear parent of node "<<print()<<endl;;

		parentEdge = NULL;
		minSlackEdges.clear();
		slack = std::numeric_limits<double>::infinity();
		inTree = false;
}


bool BiNode::isInTree() {
	return inTree;
}
void BiNode::putInTree() {
		if(inTree ) {
			cout<<"Try to put in tree a node already in the tree."<<endl;
			exit(1);
		}
		inTree = true;
}

void BiNode::removeFromTree() {
	if(inTree == false ) {
			cout<<"Try to remove from the tree a node not in the tree."<<endl;
			exit(1);
		}
		inTree = false;
}


double BiNode::getSlack() const {
		if(minSlackEdges.size() == 0 ) {
			return std::numeric_limits<double>::infinity();
		}
		return minSlackEdges.front()->getReducedCost();
}


void BiNode::addMinSlackEdge(BiEdge* edge) {
	if(minSlackEdges.size() > 0) {
		if(edge->getReducedCost() < minSlackEdges.front()->getReducedCost() - PRECISION) {
			minSlackEdges.clear();
		}
		else if ( edge->getReducedCost() > minSlackEdges.front()->getReducedCost() + PRECISION ) {
			cout<<"Error, trying to add a min slack edges which reduced cost is greater than slack value."<<endl;
			exit(1);
		}
	}
	minSlackEdges.push_back(edge);
}

bigraphEdges* BiNode::getMinSlackEdges() {
		return &minSlackEdges;
}



void BiNode::init() {
	if(edges->size() == 0) {
			cout<<"Error while initializing slacks for hungarian method"<<endl;
			cout<<"Found a node with no edge. "<<print()<<endl;
			exit(1);
	}
	
	//clear slack, min slack edges an inTree
	clear();
	potential = 0.0;
	matchingEdge = NULL;
}


list<BiNode*>::iterator BiNode::getPosition() {
		return position_listNotInTree;
}
void BiNode::setPosition(list<BiNode*>::iterator pos) {
		position_listNotInTree = pos;
}

BiEdge* BiNode::getMatchingEdge() {
	return matchingEdge;
}
void BiNode::setMatchingEdge(BiEdge * edge) {
	if(set) {
		if(edge->getLeftNode()->getIndex() != index) {
			cout<<"Error. Try to match node "<<index<<" using edge "<<edge->print()<<endl;
			exit(1);
		}
	}
	else {
		if(edge->getRightNode()->getIndex() != index) {
			cout<<"Error. Try to match node "<<index<<" using edge "<<edge->print()<<endl;
			exit(1);
		}
	}
	matchingEdge = edge;
}

string BiNode::print() const {
/*
	cout<<"printNode";
	cout.flush();
	
	*/
	std::stringstream str(stringstream::in | stringstream::out);
	
	/*
	cout<<" uh set is "<<set<<endl;
		cout.flush();
	*/
	if(set) {

		str.put('L');

		

	}
	else {

		str.put('R');


	}
	//cout<<";;;;done"<<endl;
	//cout.flush();
	str <<index;
	
	return str.str();
}

