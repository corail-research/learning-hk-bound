#include <vector>
#include <definitions.h>
#include <parse.h>


/* global parameters */
extern const char* FileIn;
extern const char* FileOut;
extern const char* FileSol;
extern const char* dotFile;

extern bool isSym;
extern bool tour_file_specified;
extern bool write_output;
extern int output_level;
extern bool prim;
extern bool prim2;
extern bool kruskal;
extern bool printDot;
extern double upperBound;
extern bool upperBoundSpecified;
extern bool heldKarp;
extern bool filter;
extern bool heldKarp2;
extern bool filter2;

extern bool forceCost;
extern bool forceCut;
extern bool forceDegree;
extern int nbForce;

extern bool filterDegree;

extern bool heavyFiltering;
extern bool noFiltering;

extern bool printReducedGraph;
extern const char* reduceGraphFileName;

extern bool useAPBound;
extern bool filterAP;

extern bool ap_bb;

/* Held-Karp parameters */
extern bool boost;

//////////////////////////////////////////////////

void usage() {
  cout << endl << "Usage: ./tsp -filein <file> [options]" << endl 
	<< "where input file must be in TSPLIB format" << endl << endl
	<< "I/O Options:" << endl 
	<< "  -fileout <file>  to specify output file to write some information " << endl
	<<endl 
	<< "  -fileSol <file>  to specify a tour file " << endl
	<< endl 
	<< "  -printDot <outfile>  print the graph in dot format in outfile " <<endl    
	<< endl
	<< "  -printReducedGraph <outfile>  print the graph obtained after a round of filtering in TSPLIB format " <<endl    
	<< endl
	<<endl << "Branch and Bound Options:" << endl 
	<< "  -useAPBound  use assignment problem bound on ATSP." <<endl    
	<< endl
	
	/*
	<< "MST Options:" << endl 
	<< "  -Prim   compute the minimum spanning tree using Prim algorithm. Print cost of the tree." << endl 
	<< endl
	<< "  -Prim2   compute the minimum spanning tree using Prim algorithm with a better data structure.  Print cost of the tree." << endl 
	<< endl
	<< "  -Kruskal   compute the minimum spanning tree using Kruskal algorithm and print its cost " << endl 
	<<endl
    
	<< "Held and Karp Options:" << endl 
	<< "  -UB val   Specify val as an upper bound to the TSP problem." << endl 
	<< endl
	<< "  -hk   Compute the Held and Karp bound using the given upper bound (with -UB)." << endl 
	<< endl
	<< "  -hk2   Compute the Held and Karp bound using the given upper bound (with -UB) using another ascent scheme." << endl 
	<< endl
	<< "  -filter  Filter edges using Held and Karp as lower bound and the given upper bound (with -UB)." << endl 
	<< "           -hk or -hk2 must be called to permit filtering"<<endl
	<< endl
	<< "  -filter2  Filter edges using Held and Karp as lower bound and the given upper bound (with -UB) using a better method." << endl 
	<< "           -hk or -hk2 must be called to permit filtering"<<endl
	<< endl
	*/
	<< "  -forceCost  Force edges in the minimum spanning tree using costs methode." << endl 
	<< endl
	<< "  -forceCut  Force edges in the minimum spanning tree using number of edges in cuts." << endl 
	<< endl
	<< "  -forceDegree  Force edges connected to a node of degree 2." << endl 
	<< endl
	<< "  -nbForce  Specify how many succesive filtering/forcing are done at each branch and bound node." << endl<<"\t If set to -1, filtering/forcing are done until fixed point is reached."<<endl
	<<"Set to 1 by default."<<endl
	<<endl
	<< "  -heavyFiltering  Filtering is done at each iteration of HK" << endl
	<< endl
	<< "  -noFiltering  Desactivate all filtering" << endl
	<< endl
      
      
	<< endl<<endl
	<< "Examples:" << endl
	<< "  ./tsp -filein tsplib/gr21.tsp" <<endl
	<< "  ./tsp -filein tsplib/gr21.tsp -Prim"<<endl
	<< "  ./tsp -filein tsplib/gr21.tsp -fileSol tsplib/gr21.opt.tour" <<endl;
  exit(0);
}

void parseArgs(int argc,char **argv) 
{	
	tour_file_specified = false;
	printDot=false;
	prim=false;
	prim2=false;
	kruskal=false;
	bool input_file_specified = false;
	write_output = false;
	output_level = 0;
	nbForce=1;

	upperBoundSpecified=false;
	heldKarp=false;
	filter=false;
	heldKarp2=false;
	filter2=false;
	forceCost=false;
	forceCut=false;
	forceDegree=false;
	heavyFiltering=false;

	filterDegree = true;
	printReducedGraph=false;
	
	useAPBound=false;
	filterAP=false;
	
	ap_bb = false;
	
	cerr<<endl<<endl;
	
	
  for (int argIndex=1; argIndex < argc; ++argIndex) {
    
    if ( !strcmp(argv[argIndex], "-filein") ) {
      argIndex++;
      FileIn = argv[argIndex];
			input_file_specified = true;
		
		cout<<"-filein "<<FileIn<<"  ";
	  cerr<<"-filein "<<FileIn<<"  ";
    }
   /*
    else if ( !strcmp(argv[argIndex], "-fileout") ) {
      argIndex++;
      FileOut = argv[argIndex];
			write_output = true;
    }
    else if ( !strcmp(argv[argIndex], "-output_level") ) {
      argIndex++;
      output_level = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-fileSol") ) {
      argIndex++;
      FileSol = argv[argIndex];
		tour_file_specified = true;
		
    }
    */
    else if ( !strcmp(argv[argIndex], "-printDot") ) {
      argIndex++;
      dotFile = argv[argIndex];
		printDot = true;
		
    }
     else if ( !strcmp(argv[argIndex], "-printReducedGraph") ) {
		argIndex++;
		reduceGraphFileName = argv[argIndex];
		printReducedGraph = true;

		cout<<"-printReducedGraph "<<reduceGraphFileName<<" ";
		cerr<<"-printReducedGraph "<<reduceGraphFileName<<" ";
	}
    
    
    
    else if ( !strcmp(argv[argIndex], "-UB") ) {
      argIndex++;
      upperBound = atof(argv[argIndex]);
		upperBoundSpecified = true; 
		 cout<<"-UB "<<upperBound<<" ";
	 	 cerr<<"-UB "<<upperBound<<" ";
		
    }
     else if ( !strcmp(argv[argIndex], "-nbForce") ) {
      argIndex++;
      nbForce = atoi(argv[argIndex]);
      if(nbForce < -1 ) {
      	cout<<"error : nb force must be positive or -1 for fixed point."<<endl;
      	exit(1);
      }
		
    }
    /*
    else if ( !strcmp(argv[argIndex], "-hk") ) {
		heldKarp = true;
		
    }
    else if ( !strcmp(argv[argIndex], "-hk2") ) {
		heldKarp2 = true;
		
    }
     else if ( !strcmp(argv[argIndex], "-filter") ) {
		filter = true;
		
    }
	else if ( !strcmp(argv[argIndex], "-filter2") ) {
		filter2 = true;
		
    }
    */
    
    else if ( !strcmp(argv[argIndex], "-forceCost") ) {
		forceCost = true;
		cout<<"-forceCost"<<" ";
		cerr<<"-forceCost"<<" ";
    }
	else if ( !strcmp(argv[argIndex], "-forceCut") ) {
		forceCut = true;
		cout<<"-forceCut"<<" ";
		cerr<<"-forceCut"<<" ";
    }
    else if ( !strcmp(argv[argIndex], "-forceDegree") ) {
		forceDegree = true;
		cout<<"-forceDegree"<<" ";
		cerr<<"-forceDegree"<<" ";
    }
    else if ( !strcmp(argv[argIndex], "-heavyFiltering") ) {
		heavyFiltering = true;
		cout<<"-heavyFiltering"<<" ";
		cerr<<"-heavyFiltering"<<" ";
    }
    else if ( !strcmp(argv[argIndex], "-noFiltering") ) {
		noFiltering=true;
		forceCost = false;
		cout<<"-noFiltering"<<" ";
		cerr<<"-noFiltering"<<" ";
    }
    
    
    
     else if ( !strcmp(argv[argIndex], "-useAPBound") ) {
		useAPBound = true;
		cout<<"-useAPBound"<<" ";
		cerr<<"-useAPBound"<<" ";
    }
	else if ( !strcmp(argv[argIndex], "-filterAP") ) {
		filterAP = true;
		useAPBound = true;
		cout<<"-useAPBound"<<" ";
		cerr<<"-useAPBound"<<" ";
		cout<<"-filterAP"<<" ";
		cerr<<"-filterAP"<<" ";
    }
	
	 else if ( !strcmp(argv[argIndex], "-AP") ) {
		ap_bb = true;
		cout<<"-AP"<<" ";
		cerr<<"-AP"<<" ";
    }
    
    /*
    else if ( !strcmp(argv[argIndex], "-Prim") ) {
		prim = true;
		
    }
    else if ( !strcmp(argv[argIndex], "-Prim2") ) {
		prim2 = true;
		
    }
    else if ( !strcmp(argv[argIndex], "-Kruskal") ) {
		kruskal = true;
		
    }
  */
    else if ( !strcmp(argv[argIndex], "-h") ) {
			usage();
    }
    else if ( !strcmp(argv[argIndex], "-help") ) {
			usage();
    }
    else {
      cout << "Unexpected option: " << argv[argIndex] << endl;
			//usage();
    }
  }
	if (!input_file_specified) {
		cout << "Please specify an input file" << endl;
		//usage();
	}
	/*
	if(filter && !heldKarp && ! heldKarp2) {
		cout << "-filter must be called along with -hk or -hk2" << endl;
		exit(1);
	}
	if(filter2 && !heldKarp && ! heldKarp2) {
		cout << "-filter2 must be called along with -hk or -hk2" << endl;
		exit(1);
	}
	*/
	cerr<<"-nbForce "<<nbForce;
	cerr<<endl;
}

//////////////////////////////////////////////////

Graphe* readTSPInstance(int& size) {
  //COULD ALSO PASS DATA STRUCTURE AS ARGUMENT
  //readTSPInstance(int&n, YOUR_DATA_STRUCTURE);

  

  ifstream inFile(FileIn, ios::in);
  if (!inFile.is_open()) {
    cout << "error while opening file " << FileIn << endl;
    exit(1);
  }
  
  char foo[255];
  char format[255];
  int i,j;
  
  char type[255];
  //read type
  while(strcmp(foo,"TYPE:") && strcmp(foo,"TYPE")) {
    inFile >> foo;
    if(inFile.eof()) {
    	cout << "End of file reached while waiting for TYPE "<<endl;
    	exit(0);
    }
  }
  
  if (!strcmp(foo,"TYPE:"))
    inFile >> type;
  else {
    inFile >> foo; // read the :
    inFile >> type;

  }
  
  //cout << "type : "<<type <<endl;
  
  //read dimension
  while(strcmp(foo,"DIMENSION:") && strcmp(foo,"DIMENSION")) {
    inFile >> foo;
    if(inFile.eof()) {
    	cout << "End of file reached while waiting for DIMENSION "<<endl;
    	exit(0);
    }
  }
  
  if (strcmp(foo,"DIMENSION:") && strcmp(foo,"DIMENSION")) {
    cout << "could not read file!" << endl;
    exit(0);
  }
  
  if (!strcmp(foo,"DIMENSION:"))
    inFile >> size;
  else {
    inFile >> foo; // read the :
    inFile >> size;

  }
  
	Graphe * graphe;
	
	if(!strcmp(type,"TSP")) {
		graphe = new SymTSP(size, false);
		isSym=true;
	}
	else if(!strcmp(type,"ATSP")) {
		graphe = new ATSP(size);
		isSym=false;
	}
	else {
		cout <<"type "<<type<<" not supported"<<endl;
		exit(1);
	}

 
 // cout << "DIMENSION: " << size << endl;
  //cout << size <<"\t";
	while(strcmp(foo,"EDGE_WEIGHT_FORMAT:") && strcmp(foo,"EDGE_WEIGHT_FORMAT"))  {
		inFile >> foo;
		if(inFile.eof()) {
    		cout << "End of file reached while waiting for EDGE_WEIGHT_FORMAT "<<endl;
    		exit(0);
   		}
	}
	if ( strcmp(foo,"EDGE_WEIGHT_FORMAT:") && strcmp(foo,"EDGE_WEIGHT_FORMAT")) {
	  cout << "could not read file!" << endl;
	  exit(0);
	}
	if (!strcmp(foo,"EDGE_WEIGHT_FORMAT:")) 
		inFile >> format;
	else {
		inFile >> foo; // read :
		inFile >> format;
	}
		
	//cout << "format = " << format << endl;
	if (! strcmp(format,"LOWER_DIAG_ROW")) {

		while(strcmp(foo,"EDGE_WEIGHT_SECTION")) { 
			inFile >> foo;
			if(inFile.eof()) {
				cout << "End of file reached while waiting for EDGE_WEIGHT_SECTION "<<endl;
				exit(0);
	   		}	
		}
		if (strcmp(foo,"EDGE_WEIGHT_SECTION")) {
			cout << "could not read file!" << endl;
			exit(0);
		}
		for(i = 0; i < size; i++){
			for(j = 0; j <= i; j++) {
				double foo;
				inFile >> foo;
//			cout << "i=" << i << ",j=" << j << " foo" << foo << endl;
				if ((foo < NO_EDGE) && (i!=j)) {
//				cout << "add to edges" << endl;
				  //ADD AN EDGE HERE
					graphe->addEdge(i,j,foo);
				}
			}
		}
	}
	else if ((! strcmp(format,"FULL_DIST")) || (! strcmp(format,"FULL_MATRIX"))) {
		while(strcmp(foo,"EDGE_WEIGHT_SECTION")) { 
			inFile >> foo ;
			
			if(inFile.eof()) {
				cout << "End of file reached while waiting for EDGE_WEIGHT_SECTION "<<endl;
				exit(0);
	   		}	
		}
		if (strcmp(foo,"EDGE_WEIGHT_SECTION")) {
			cout << "could not read file!" << endl;
			exit(0);
		}
		for(i = 0; i < size ; i++) {
			for(j = 0; j < size; j++) {
				double foo;
				inFile >> foo;
				if ((foo < NO_EDGE) && (i!=j)) {
					//ADD AN EDGE HERE
					if(!isSym) {
						graphe->addEdge(i,j,foo);
					}
					else if(j<i){
						graphe->addEdge(i,j,foo);
					}
				}
			}
		}
	}
	else if (! strcmp(format,"UPPER_ROW")) { // first line contains n-1 entries, second line n-2, etc... i==j means inf
		while(strcmp(foo,"EDGE_WEIGHT_SECTION")) {
			inFile >> foo ;
			if(inFile.eof()) {
				cout << "End of file reached while waiting for EDGE_WEIGHT_SECTION "<<endl;
				exit(0);
	   		}	
		}
		if (strcmp(foo,"EDGE_WEIGHT_SECTION")) {
			cout << "could not read file!" << endl;
			exit(0);
		}
		for(i = 0; i < size-1; i++) {
			for(j = i+1; j < size; j++) {
				double foo;
				inFile >> foo;
				if ((foo < NO_EDGE) && (i!=j)) {
				  //ADD AN EDGE HERE
				  //cout<< i << j<< endl;
					graphe->addEdge(i,j,foo);
				}
			}
		}
	}
	else if (! strcmp(format,"UPPER_DIAG_ROW")) {
		cout << "UPPER_DIAG_ROW not yet supported" << endl;
		exit(0);
	}
	else {
		cout <<"format "<<format <<"not recognized"<<endl;
		exit(0);
	}
	

	inFile.close();
	return graphe;
}




Tour* readTourInstance() {
	if(!tour_file_specified) {
		cout << "no tour file specified " << endl;
		cout << "use -fileSol fileName to specify one" << endl ;
		exit(1);
	}


	ifstream solFile(FileSol, ios::in);
	if (!solFile.is_open()) {
	cout << "error while opening file " << FileIn << endl;
	exit(1);
	}

	char foo[255];

	//read the dimension

	while(strcmp(foo,"DIMENSION:") && strcmp(foo,"DIMENSION")) {
	solFile >> foo;
	
	}
	if (strcmp(foo,"DIMENSION:") && strcmp(foo,"DIMENSION")) {
	cout << "could not read file!" << endl;
	exit(0);
	}

	int size;

	if (!strcmp(foo,"DIMENSION:"))
		solFile >> size;
	else {
		solFile >> foo; // read the :
		solFile >> size;

	}
  	cout << "size of tour is : "<<size << endl;
	//Tour * tour = new Tour(size);
  
  
	while(strcmp(foo,"TOUR_SECTION")) solFile >> foo;
	
	if (strcmp(foo,"TOUR_SECTION")) {
		cout << "could not read file!" << endl;
		exit(0);
	}
	
	vector<int> tour;
	tour.reserve(size);
	for( int i=0; i<size; i++) {
		int node;
		solFile>>node;
		tour.push_back(node-1);
	}
	
	int lastNode;
	solFile>>lastNode;
	
	if(lastNode != -1) {
		cout << "Tour section must end with a -1" <<endl;
		exit(1);
	}
  
	solFile.close();
	
	return new Tour(tour);
}

Graphe* readInputFile() {

	

  enum distanceType { unknownDistance, weights, coords };
  enum weightFormat { unknownWeight, full, lower, upper, lowerdiag, upperdiag };
  enum edgeWeightType {unknownType, explici ,euc_2d, max_2d,  man_2d, ceil_2d, geo, att };

  ifstream inFile(FileIn, ios::in);
  if (!inFile.is_open()) {
    cout << "error while opening file " << FileIn << endl;
    exit(1);
  }
	int NNodes;  
  char foo[255];
  char type[255];
  bool foundDimension = false;
  bool foundType=false;
  
  distanceType DType = unknownDistance;
  weightFormat WFormat = unknownWeight;
  edgeWeightType WType= unknownType;

  while (!inFile.eof()) {
    inFile >> foo;
    if ( (!strcmp(foo, "DIMENSION")) || (!strcmp(foo, "dimension")) ) {
      inFile >> foo; // read the additional " : "
      inFile >> NNodes;
      foundDimension = true;
    }
    else if ( (!strcmp(foo, "DIMENSION:")) || (!strcmp(foo, "dimension:")) ) {
      inFile >> NNodes;
      foundDimension = true;
    }
    else if ( (!strcmp(foo, "NODE_COORD_SECTION")) || (!strcmp(foo, "node_coord_section")) ) {
      DType = coords;
    }
    else if ( (!strcmp(foo, "EDGE_WEIGHT_SECTION")) || (!strcmp(foo, "edge_weight_section")) ) {
      DType = weights;
    }
    else if ( (!strcmp(foo, "LOWER_ROW")) || (!strcmp(foo, "lower_row")) ) {
      WFormat = lower;
    }
    else if ( (!strcmp(foo, "UPPER_ROW")) || (!strcmp(foo, "upper_row")) ) {
      WFormat = upper;
    }
    else if ( (!strcmp(foo, "LOWER_DIAG_ROW")) || (!strcmp(foo, "lower_diag_row")) ) {
      WFormat = lowerdiag;
    }
    else if ( (!strcmp(foo, "UPPER_DIAG_ROW")) || (!strcmp(foo, "upper_diag_row")) ) {
      WFormat = upperdiag;
    }
    else if ( (!strcmp(foo, "FULL_MATRIX")) || (!strcmp(foo, "full_matrix")) ) {
      WFormat = full;
    }
    
     else if ( (!strcmp(foo, "EXPLICIT")) || (!strcmp(foo, "explicit")) ) {
      WType = explici; 
    }
     else if ( (!strcmp(foo, "EUC_2D")) || (!strcmp(foo, "euc_2d")) ) {
      WType = euc_2d; 
    }
	else if ( (!strcmp(foo, "MAX_2D")) || (!strcmp(foo, "max_2d")) ) {
      WType = max_2d; 
    }
	else if ( (!strcmp(foo, "MAN_2D")) || (!strcmp(foo, "man_2d")) ) {
      WType = man_2d; 
    }
	else if ( (!strcmp(foo, "CEIL_2D")) || (!strcmp(foo, "ceil_2d")) ) {
      WType = ceil_2d; 
    }
    else if ( (!strcmp(foo, "GEO")) || (!strcmp(foo, "geo")) ) {
      WType = geo; 
    }
    else if ( (!strcmp(foo, "ATT")) || (!strcmp(foo, "att")) ) {
      WType = att; 
    }

    
    
    
    else if (!strcmp(foo,"TYPE:") ) {
    	inFile >> type;
    	foundType=true;
    }
	else if (!strcmp(foo,"TYPE") ){
		inFile >> foo; // read the :
		inFile >> type;
		foundType=true;
	}

    if (foundDimension && !(DType == unknownDistance) && !(WFormat == unknownWeight) && WType != unknownType  &&foundType) {
      break;
    } 
  }

  inFile.close();

  if (!foundDimension) {
    cout << "error: could not determine dimension while reading file" << endl;
    exit(1);
  }
  if (DType == unknownDistance) {
    cout << "error: could not determine distance type while reading file" << endl;
    exit(1);
  }
  if ((DType == weights) && (WFormat == unknownWeight)) {
    cout << "error: could not determine format of the edge weights while reading file" << endl;
    exit(1);
  }
  if(WType == unknownType) {
  	cout<<"error : edge_weight_type not supported"<<endl;
  	exit(1);
  }
  if (!foundType) {
  	cout<<"error : no type found in file (stsp or atps)."<<endl;
  	exit(1);
  }
  
 
  Graphe * graphe;
	if(!strcmp(type,"TSP")) {
		graphe = new SymTSP(NNodes, false);
		isSym=true;
	}
	else if(!strcmp(type,"ATSP")) {
		graphe = new ATSP(NNodes);
		isSym=false;
	}
	else {
		cout <<"type "<<type<<" not supported"<<endl;
		exit(1);
	}

/*
  Distances = new double*[NNodes];
  for (int i=0; i<NNodes; i++) {
    Distances[i] = new double[NNodes];
    for (int j=0; j<NNodes; j++) {
      Distances[i][j] = 0;
    }
  }
*/
  double** Coordinates;
  if (DType == coords) {
    Coordinates = new double*[NNodes];
    for (int i=0; i<NNodes; i++) {
      Coordinates[i] = new double[2]; // [x,y] coordinates
      for (int j=0; j<2; j++) {
	Coordinates[i][j] = 0;
      }
    }    
  }

  inFile.open(FileIn, ios::in);

  if (DType == weights) {
    if (WFormat == full) {
      inFile >> foo;
      while ( (strcmp(foo, "EDGE_WEIGHT_SECTION")) && (strcmp(foo, "edge_weight_section")) ) {
	inFile >> foo;
	if (inFile.eof()) {
	  cout << "error while reading distance matrix" << endl;
	  exit(1);
	}
      }
      for (int i=0; i<NNodes; i++) {
	for (int j=0; j<NNodes; j++) {
	
	
//	  inFile >> Distances[i][j];
		double foo;
		inFile >> foo;
		if( j > i || (j <i && !isSym) ) {
			if( foo < NO_EDGE) {
				graphe->addEdge(i,j,foo);
			}
		}


	  if (inFile.eof() && (i!=NNodes-1) && (j!=NNodes-1)) {
	    cout << "error while reading distance matrix:" << endl
		 << "end of file reached while reading entry (" << i << "," << j << ")" << endl;
	    exit(1);
	  }
	}
      }
    }
    else if ((WFormat == upper) || (WFormat == upperdiag)) {
      inFile >> foo;
      while ( (strcmp(foo, "EDGE_WEIGHT_SECTION")) && (strcmp(foo, "edge_weight_section")) ) {
	inFile >> foo;
	if (inFile.eof()) {
	  cout << "error while reading distance matrix" << endl;
	  exit(1);
	}
      }
	for (int i=0; i<NNodes; i++) {
		for (int j=i; j<NNodes; j++) {
			if (i==j) {
				if (WFormat == upperdiag) {
				  double foo;
				  inFile >> foo;
				 // graphe->addEdge(i,j,foo);
				}
				else {
				  //Distances[i][j] = 0;
				}
			}
			else {
			 double foo;
			  inFile >> foo;
			  if( foo < NO_EDGE) {
				  graphe->addEdge(i,j,foo);
				  if( !isSym) {
			  			graphe->addEdge(j,i,foo);
				  }
			  }
			//Distances[j][i] = Distances[i][j];
				if (inFile.eof() && (i!=NNodes-2) && (j!=NNodes-2)) {
				  cout << "error while reading distance matrix:" << endl
				   << "end of file reached while reading entry (" << i << "," << j << ")" << endl;
				  exit(1);
				}
			}
		}
	}
    }
    else if ((WFormat == lower) || (WFormat == lowerdiag))  {
      inFile >> foo;
      while ( (strcmp(foo, "EDGE_WEIGHT_SECTION")) && (strcmp(foo, "edge_weight_section")) ) {
	inFile >> foo;
	if (inFile.eof()) {
	  cout << "error while reading distance matrix" << endl;
	  exit(1);
	}
      }
      for (int j=0; j<NNodes; j++) {
	for (int i=0; i<=j; i++) {
	  if (i==j) {

	  	
	    if (WFormat == lowerdiag) {
			double foo;
			inFile >> foo;
			//graphe->addEdge(i,j,foo);
	    }
	    else {
	    //  Distances[i][j] = 0;
	    }
	    
	  }
	  else {
		double foo;
		inFile >> foo;
		if( foo < NO_EDGE) {
			graphe->addEdge(i,j,foo);
			if( !isSym) {
				graphe->addEdge(j,i,foo);
			}
		}
	   // Distances[j][i] = Distances[i][j];
	    if (inFile.eof() && (i!=NNodes-2) && (j!=NNodes-2)) {
	      cout << "error while reading full distance matrix:" << endl
		   << "end of file reached while reading entry (" << i << "," << j << ")" << endl;
	      exit(1);
	    }
	  }
	}
      }
    }
  }
  else if (DType == coords) {
    inFile >> foo;
    while ( (strcmp(foo, "NODE_COORD_SECTION")) && (strcmp(foo, "node_coord_section")) ) {
      inFile >> foo;
      if (inFile.eof()) {
	cout << "error while reading coordinate matrix" << endl;
	exit(1);
      }
    }
    for (int i=0; i<NNodes; i++) {
      int tmp;
      inFile >> tmp;
      if (tmp != i+1) {
	cout << "error while reading coordinates of node " << i+1 << endl;
	exit(1);
      }
      inFile >> Coordinates[i][0];
      inFile >> Coordinates[i][1];
    //  cout<<"node "<<i<<" "<<Coordinates[i][0]<<" "<<Coordinates[i][1]<<endl;
      if (inFile.eof() && (i!=NNodes-1)) {
	cout << "error while reading distance matrix:" << endl
	     << "end of file reached while reading entry " << i << endl;
	exit(1);
      }
    }
    for (int i=0; i<NNodes; i++) {
      for (int j=i+1; j<NNodes; j++) {
      	double foo;
      	if(WType==euc_2d) {
      		double xd= Coordinates[i][0] - Coordinates[j][0];
      		double yd= Coordinates[i][1] - Coordinates[j][1];
      		foo = round(sqrt(xd*xd+yd*yd));
      	}
      	else if (WType == man_2d) {
      		double xd= fabs(Coordinates[i][0] - Coordinates[j][0]);
      		double yd= fabs(Coordinates[i][1] - Coordinates[j][1]);
      		foo = round(xd+yd);
      	}
      	else if (WType == max_2d) {
      		double xd= round( fabs(Coordinates[i][0] - Coordinates[j][0]) );
      		double yd= round( fabs(Coordinates[i][1] - Coordinates[j][1]) );
      		if(yd > xd) {
      			foo = yd;
      		}
      		else {
      			foo=xd;
      		}	
      	}
      	else if (WType == geo) {
      		const double PI = 3.141592;
      		
      		int deg= floor(Coordinates[i][0]);
      		double min = Coordinates[i][0] - deg;
      		double latitude_i= PI * (deg + 5.0* min / 3.0) / 180.0;
      		
      		deg = floor (Coordinates[i][1] );
      		min = Coordinates[i][1] - deg;
      		double longitude_i= PI * (deg + 5.0 * min / 3.0 ) / 180.0;
      		
      		deg = floor (Coordinates[j][0] );
      		min = Coordinates[j][0] - deg;
      		double latitude_j= PI * (deg + 5.0 * min / 3.0 ) / 180.0;
      		
      		deg = floor (Coordinates[j][1] );
      		min = Coordinates[j][1] - deg;
      		double longitude_j= PI * (deg + 5.0 * min / 3.0 ) / 180.0;
      		
      		const double RRR = 6378.388;
      		
      		double q1 = cos( longitude_i - longitude_j);
      		double q2 = cos( latitude_i - latitude_j);
      		double q3 = cos( latitude_i + latitude_j);
      		
      		foo = (int) ( RRR * acos( 0.5 * ((1.0+q1)*q2 - (1.0-q1)*q3) ) +1.0 );
      	}
      	else if (WType == att) {
      		double xd= Coordinates[i][0] - Coordinates[j][0];
      		double yd= Coordinates[i][1] - Coordinates[j][1];
      		double r = sqrt((xd*xd+yd*yd)/10.0);
      		double t = round(r);
      		if(t<r) {
      			foo = t+1;
  			}
  			else {
  				foo=t;
  			}
      	}
      	else if(WType==ceil_2d) {
      		double xd= Coordinates[i][0] - Coordinates[j][0];
      		double yd= Coordinates[i][1] - Coordinates[j][1];
      		foo = ceil(sqrt(xd*xd+yd*yd));
      	}
      	else {
      		cout<<"Error in parsing, no edge weight type specified"<<endl;
      		exit(1);
      	}

	       graphe->addEdge(i,j,foo);
	       
	       if(!isSym) {
	       		graphe->addEdge(j,i,foo);
	       }
//	Distances[j][i] = Distances[i][j];
      }
    }
  }
  else {
    cout << "error while reading data: unknown distance type" << endl;
    exit(1);
  }

 
  return graphe;
}

