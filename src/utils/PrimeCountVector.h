#ifndef PRIMECOUNTVECTOR_H
#define PRIMECOUNTVECTOR_H

#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>
#include <algorithm>

namespace primecount {

template <typename T>
class PrimeCountVector : public std::vector<T>,
                         public PrimeSieveCallback<uint64_t> {
public:
  void callback(uint64_t prime)
  {
    this->push_back( static_cast<T>(prime) );
  }
  /// Count primes using binary search O(log(n))
  T pi2(T n) const
  {
    return static_cast<T>(
        std::upper_bound(this->begin(), this->end(), n) - this->begin()
    );
  }
                           T pi(T schluessel) const
                           {
                             int links = 0;                  // linke Teilfeldbegrenzung
                             int rechts = this->size() - 1;  // rechte Teilfeldbegrenzung
                             int versch;                     // Anzahl verschiedener Elemente
                             int pos;                        // aktuelle Teilungsposition

                             // solange der Schlüssel im Bereich liegt (andernfalls ist das gesuchte
                             // Element nicht vorhanden)
                             while(this->operator[](links) <= schluessel && this->operator[](rechts) > schluessel){
                               // Aktualisierung der Anzahl der verschiedenen Elemente
                               versch = this->operator[](rechts) - this->operator[](links);
                               
                               // Berechnung der neuen interpolierten Teilungsposition
                               pos = links + (int)(((double)rechts - links) * (schluessel - this->operator[](links))
                                                   / versch);
                               
                               if(this->operator[](pos) <= schluessel)            // rechtes Teilintervall
                                 links = pos + 1;                       // daten[pos] bereits überprüft
                               else if(this->operator[](pos) > schluessel)      // linkes Teilintervall
                                 rechts = pos - 1;                      // daten[pos] bereits überprüft
                               else                                     // Element gefunden
                                 return pos;                            // Position zurückgeben
                             }
                             
                             if (this->operator[](pos) <= schluessel)
                               return pos + 1;

                             return pos;                                 // Element nicht gefunden
                           }
};
  
} // end namespace

#endif
