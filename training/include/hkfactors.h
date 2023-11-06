//////////////////////////////////////////////////
#define DBG(x)       // use this to not print debug info
//#define DBG(x) x     // use this to print debug info


#include <definitions.h>
#include <parse.h>
#include <hk.h>
#include <bb.h>
#include <bigraphe.h>
#include <hungarianMethod.h>


extern "C"{

float *hkfactors(int argc, char* argv[]);

void hkfact_del(float *bound_factors);

}
