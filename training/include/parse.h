#include <graphe.h>
#include <tour.h>
#include <stsp.h>
#include <atsp.h>
#include <math.h>

/**
 *
 * \brief Parse command arguments.
 */
void parseArgs(int argc,char **argv);
/**
 *
 * \brief Read a TSP instance from the given input file and return the corresponding Graphe object.
 * Return a SymTSP or a ATSP object according to the field TYPE in the given input file.
 *
 */
Graphe* readTSPInstance(int& size);
/**
 *
 * \brief Read a tour from the given input file and return the correponding Tour object.
 */
Tour* readTourInstance();

Graphe* readInputFile();

