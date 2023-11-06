#include <node.h> 

Node::Node() :
	index(NULL),
	nbForcedVoisins(0)
{}

Node::Node(int i) :
	index(i),
	nbForcedVoisins(0)
{
	return; 
}

nodeEdges::iterator Node::addEdge( Edge * e) {
	const int source = e->getSource()->getIndex();
	const int  dest = e->getDest()->getIndex();
	
	// Check if e can be added to this node
	if ( dest!=getIndex() && source!=getIndex() ) {
		cout << "Bad edge";
		exit(1);
	}
	if(source==dest) {
		cout << "Bad edge";
		exit(1);
	}
	
	/*
	// Remember the new neighboor in voisins.
	if(dest==getIndex()) {
		voisins.insert(dest);
	}	
	else {
		voisins.insert(source);
	}
	*/
	//Remember e in edges.
	edges.push_back(e);
	
	return --(edges.end());
		
}

void Node::removeEdge(nodeEdges::iterator position){
	edges.erase(position);
}


void Node::addForcedVoisin() {
	nbForcedVoisins++;
}
void Node::removeForcedVoisin() {
	nbForcedVoisins--;
	if(nbForcedVoisins <0) {
		cout<<"Nb forced voisins =0"<<endl;
		exit(1);
	}
}

int Node::getNbForcedVoisins() {
	return nbForcedVoisins;
}

int Node::getIndex() const {
	int i=index;
	return i;
}

/*
const nodeVoisins* Node::getVoisins() {
 	return &voisins;   
}
*/
nodeEdges* Node::getEdges() {
	return &edges;
}

bool Node::operator<(const Node & b) const {
	return this->index < b.index;
}

bool  Node::operator==(const Node & b) const {
	return this->index == b.index;
}
bool  Node::operator!=(const Node & b) const {
	return this->index != b.index;
}
