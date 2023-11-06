#ifndef TOUR_H
#define TOUR_H
#include <vector>
#include <definitions.h>
/*!
 * \class Tour
 * \brief A tour for a TSP instance.
 * \author Pascal Benchimol
 */


class Tour {
private :
const std::vector<int> tour; /*!< Sequence of nodes in the tour*/
public :
Tour( std::vector<int> list);

/*!
*  \brief Return the number of nodes in this tour.
*/
int getSize() const;
/*!
*  \brief Return the sequence of nodes' indexes of the tour.
*/
const std::vector<int> * getTour() const;
void print() const;
};
#endif
