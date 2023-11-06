#include <heap.h>

int Heap::getParent(int i) {
	if( i<=0)
	return 0;
	
	if(i> tableau.size()) {
		std::cout <<"Can't access parent of element "<<i<<" in a heap of size "<<tableau.size();
	}

	return floor((i-1)/2);
}


void Heap::monterDsHeap(Pointeur * p) {

	bool bienPlace=false;
	
	int i=p->tabIndex;
	
	if( i<0 || i>= tableau.size()) {
		cout <<"Error in heap"<<endl;
		cout << "asked to move up a pointeur which index is "<<i<<" in an array of size "<<tableau.size()<<endl;
		cout << "node index is "<< p->nodeRef->getIndex()<<endl;
		exit(1);
	}
	

	while(i >0 && !bienPlace) {
		int pIndex = floor((i-1)/2);
		Pointeur * parent = *(tableau.begin()+pIndex); 
		if(*parent < *p) {
			bienPlace=true;
		}
		else {
			tableau[i]=parent;
			parent->tabIndex=i;
			i=floor((i-1)/2);
			
		}
	}
	tableau[i]=p;
	p->tabIndex=i;

}

void Heap::descendreDsHeap(Pointeur* p) {


	
//	Pointeur * toInsert=tableau.front();
	int i=p->tabIndex;
	
	if( i<0 || i>= tableau.size()) {
		cout << "asked to move down a pointeur which index is "<<i<<" in an array of size "<<tableau.size()<<endl;
		exit(1);
	}

	
	
	bool bienPlace=false;
	bool noFils=false;

	while(!noFils && !bienPlace){
	
		if(2*i+2 < tableau.size() ) {
			Pointeur * filsG = *(tableau.begin()+2*i+1);
			Pointeur * filsD = *(tableau.begin()+2*i+2);
		
			//Si parent <= fils, on arete
			if( *p < *filsG && *p < *filsD ) {
				bienPlace=true;
			}
			//Sinon on remonte le fils le plus petit
			else if ( *filsG < *filsD ) {
				tableau[i]=filsG;
				filsG->tabIndex=i;
				i=2*i+1;
			}
			//Sinon filsD <= filsG
			else {
				tableau[i]=filsD;
				filsD->tabIndex=i;
				i=2*i+2;
			}
		}
		
		//Si pas de fils droit
		else if (2*i+1 < tableau.size() ) {
			Pointeur * filsG = *(tableau.begin()+2*i+1);
			//Si parent <= fils, on arete
			if( *p < *filsG) {
				bienPlace=true;
			}
			//Sinon on remonte le fils
			else {
				tableau[i]=filsG;
				filsG->tabIndex=i;
				i=2*i+1;
			}
		
		}
		
		//Si pas de fils
		else {
			noFils=true;
		}
	}
	
	tableau[i]=p;
	p->tabIndex=i;
}

void Heap::make_heap() {

	//a un fils si 2i+1 <size
	for(int i= tableau.size()/2; i>=0; --i) {
		descendreDsHeap(tableau[i]);
	}

}


Heap::Heap(std::vector<Pointeur*> & tab) :
	tableau(tab)
{	
	/*
	tableau.reserve(map.size());
	int i=0;
	for(std::map<int,Node>::iterator it = map.begin(); it!=map.end(); ++it) {
		Node* nref = &(it->second);
		Pointeur p = Pointeur(nref,i,std::numeric_limits<double>::infinity());
		tableau.push_back(&p);
	}
	*/
}

Heap::~Heap() {

}

std::vector< Pointeur*> * Heap::getTab() {
	return &tableau;
}


Pointeur * Heap::pop_heap() {


	Pointeur * result = tableau.front();
	
	tableau.front()=tableau.back();
	tableau.front()->tabIndex=0;
	tableau.pop_back();
	
	//descendre l'objet dans le heap
	if(tableau.size()>0) {
		descendreDsHeap(tableau.front());	
	}
	return result;
}



void Heap::push_heap( Pointeur * p) {

	p->tabIndex=tableau.size();	
	tableau.push_back(p);
	monterDsHeap(p);
	
}

void Heap::print() {
	for(std::vector<Pointeur*>::iterator it = tableau.begin(); it!= tableau.end(); ++it ) {
		cout << (*it)->nodeRef->getIndex()<<" ";
	}
	cout << endl;
}




