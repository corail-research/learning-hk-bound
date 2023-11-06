#include <atsp.h>

ATSP::ATSP(int n) : Graphe(n)
{}

SymTSP * ATSP::getTSP() {
	SymTSP * newGraph = new SymTSP(getSize()*2, true);
	//in node : i
	//out node :size+i
	bool canContainTour = true;
	for(int i=0; i<getSize(); i++) {
		Edge * e=newGraph->addEdge(i,i+getSize(),0.0);
		
		e->force(&canContainTour);
		
		if(canContainTour == false) {
			cout<<"Error while creating STSP from ATSP.";
			cout<<"There is 3 or more mandatory edges on one node"<<endl;
			exit(1);
		}
	}

	for(graphEdges::iterator it = edges.begin(); it!=edges.end(); ++it) {
		int newSource = (*it)->getSource()->getIndex() +getSize();
		int newDest= (*it)->getDest()->getIndex();
		newGraph->addEdge(newSource,newDest, (*it)->getWeight());
	}

	return newGraph;

}


void ATSP::printDot(string nomFichier) const {
// Ouvre un fichier en écriture
    ofstream fb(nomFichier.c_str());
    // Teste si le fichier est ouvert :
    if (fb.is_open())
    {
    	fb <<"digraph g{"<<endl;
    
      
    	
		for( graphEdges::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
			//Edge ed = *it;
			fb << (*it)->getSource()->getIndex() <<" -> " << (*it)->getDest()->getIndex()<<
				"[label = \" "<< (*it)->getWeight()<<"\" ]" <<	endl;
			// fprintf(&fb,"%d -- %d;\n", ed.getSource()->getIndex() ,ed.getDest()->getIndex());
			// string l1=l;
			//fb.sputn(l1.data(), l1.size());
		}
		fb << "}"<<endl;


        // Ferme le fichier :
        fb.close();
    }


}
