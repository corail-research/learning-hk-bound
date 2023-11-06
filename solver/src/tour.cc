#include <tour.h>



Tour::Tour(std::vector<int> vec) :
	tour(vec)
{
}

int Tour::getSize() const {
	int s=tour.size();
	return s;
}
const std::vector<int>* Tour::getTour() const {
	return &tour;
}

void Tour::print() const {
	cout << "dimension : " << tour.size() << endl;
	int i=0;
	for(std::vector<int>::const_iterator it = tour.begin(); it!=tour.end(); ++it) {
		cout << i<< " "<< *it << endl;
		i++;
		cout.flush();
	}
}

