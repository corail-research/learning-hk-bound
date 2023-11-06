#ifndef ATSP_H
#define ATSP_H
#include <graphe.h>
#include <stsp.h>

/*!
 * \class ATSP
 * \brief An Asymmetric Travelling Salesman Problem object.
 * \author Pascal Benchimol
 */

class ATSP : public Graphe {

public :
 /*!
     *  \brief Constructeur
     *  \param n : number of nodes.
     */
ATSP(int n);

 /*!
     *  \brief Build a symmetric TSP instance from this asymmetric instance using duplication of nodes.
     *  \return A SymTSP object.
     */

SymTSP * getTSP();

void printDot(string nomFichier) const;

};
#endif
