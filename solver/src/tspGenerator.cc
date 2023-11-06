


//////////////////////////////////////////////////
#define DBG(x)       // use this to not print debug info
//#define DBG(x) x     // use this to print debug info
#include <iostream>
#include <stdlib.h>


#include <sys/time.h>



/* global parameters */
const char* fileName;
bool	out_file_specified;

int size;
bool size_specified;

int main(int argc, char* argv[]){



	if(!out_file_specified) {
		cout<<"Must specify an outfile with -outfile."<<endl;
	}
	if(!size_specified) {
		cout<<"Must specify a size with -size."<<endl;
	}
 

	// Ouvre un fichier en écriture
	ofstream fb(fileName.c_str());
	// Teste si le fichier est ouvert :
	if (fb.is_open())
	{
		fb<<"TYPE: TSP"<<endl
			<<"DIMENSION: "<<size<<endl
			<<"EDGE_WEIGHT_FORMAT: UPPER_ROW"<<endl
			<<"EDGE_WEIGHT_TYPE: EXPLICIT"<<endl
			<<"EDGE_WEIGHT_FORMAT: UPPER_ROW"<<endl
			<<"EDGE_WEIGHT_SECTION"<<endl;
			
		srand(get_seed());
		for(int i=0; i< size; i++) {
			for(int j=i+1; j < size; j++) {
				int w = (rand() % 1000);
				fb<< w << "  ";
			}
		
			fb<<endl;
		}
		
		
		 fb.close();
	}
	
	

}

unsigned long get_seed(void) {
   struct timeval tv;
  struct timezone tzp;
  gettimeofday(&tv,&tzp);
  return (( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec;

}


void parseArgs(int argc,char **argv) {
	out_file_specified=false;
	size_specified=false;
	for (int argIndex=1; argIndex < argc; ++argIndex) {

		if ( !strcmp(argv[argIndex], "-outfile") ) {
			argIndex++;
			fileName = argv[argIndex];
			out_file_specified = true;

		}
		else if ( !strcmp(argv[argIndex], "-size")  ) 	 {
			argIndex++;
			size =  atol(argv[argIndex]);
			size_specified = true;
		}
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

}


void usage() {
  cout << endl << "Usage: ./tspGenerator -outfile <file> -size <size>" << endl;
  exit(0);
}
