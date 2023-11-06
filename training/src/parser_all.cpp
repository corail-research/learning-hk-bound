#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

bool output = false; // print distance matrix

/***************************/

void usage() {
  cout << endl
       << "USAGE: parser -filein <filename> [options]" << endl 
       << "where " << endl
       << "  <filename>      is the name of the file in TSPLIB format" << endl
       << "options are " << endl
       << "  -output         to print the distance matrix" << endl
       << "Example: ./parser -filein att48.tsp -output" << endl;
  exit(1);
}


void parseArgs(int argc, char **argv, char* &fileName) {

  bool fileNameSpecified = false;

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-filein") ) {
      argIndex++;
      fileName = argv[argIndex];
      fileNameSpecified = true;
    }
    else if ( !strcmp(argv[argIndex], "-h") ) {
      usage();
    }
    else if ( !strcmp(argv[argIndex], "-output") ) {
      output = true;
    }
    else {
      cout << "Unexpected option: " << argv[argIndex] << endl;
      usage();
    }
  }
  if (!fileNameSpecified) {
    cout << "please specify an input file" << endl;
    usage();
  }
}

/***************************/

Graphe* readInputFile(const char* fileName) {

	

  enum distanceType { unknownDistance, weights, coords };
  enum weightFormat { unknownWeight, full, lower, upper, lowerdiag, upperdiag };

  ifstream inFile(fileName, ios::in);

	int NNodes;  
  char foo[255];
  char type[255];
  bool foundDimension = false;
  bool foundType=false;
  
  distanceType DType = unknownDistance;
  weightFormat WFormat = unknownWeight;

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
    else if (!strcmp(foo,"TYPE:") ) {
    	inFile >> type;
    	foundType=true;
    }
	else if (!strcmp(foo,"TYPE:") ){
		inFile >> foo; // read the :
		inFile >> type;
		foundType=true;
	}

    if (foundDimension && !(DType == unknownDistance) && !(WFormat == unknownWeight) && foundType) {
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
  
  if (!foundType) {
  	cout<<"error : no type found in file (stsp or atps)."<<endl;
  	exit(1);
  }
  
 
  Graphe * graphe;
	if(!strcmp(type,"TSP")) {
		graphe = new SymTSP(NNodes);
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

  inFile.open(fileName, ios::in);

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
		graphe->addEdge(i,j,foo);

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
				  graphe->addEdge(i,j,foo);
				}
				else {
				  //Distances[i][j] = 0;
				}
			}
			else {
			 double foo;
			  inFile >> foo;
			  graphe->addEdge(i,j,foo);
			  if( !isSym) {
		  			graphe->addEdge(j,i,foo);
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
			graphe->addEdge(i,j,foo);
	    }
	    else {
	    //  Distances[i][j] = 0;
	    }
	  }
	  else {
		double foo;
		inFile >> foo;
		graphe->addEdge(i,j,foo);
		if( !isSym) {
			graphe->addEdge(j,i,foo);
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
      if (inFile.eof() && (i!=NNodes-1)) {
	cout << "error while reading distance matrix:" << endl
	     << "end of file reached while reading entry " << i << endl;
	exit(1);
      }
    }
    for (int i=0; i<NNodes; i++) {
      for (int j=0; j<NNodes; j++) {
	double foo = sqrt( (Coordinates[i][0] - Coordinates[j][0])*(Coordinates[i][0] - Coordinates[j][0]) +
			       (Coordinates[i][1] - Coordinates[j][1])*(Coordinates[i][1] - Coordinates[j][1]) );
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

/******************************/

int main(int argc, char** argv) {

  char* fileName;
  int NNodes;
  double** Distances;

  parseArgs(argc, argv, fileName);

  cout << "reading file " << fileName << endl;
  
  bool succeed = readInputFile(fileName, NNodes, Distances);

  if (output) {
    cout << "Distances:" << endl;
    for (int i=0; i<NNodes; i++) {
      for (int j=0; j<NNodes; j++) {
	cout << Distances[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
  
  return 0;
}
