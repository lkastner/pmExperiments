/* Copyright (c) 1997-2014
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#ifndef POLYMAKE_PREPROCESS_H
#define POLYMAKE_PREPROCESS_H

#include "polymake/client.h"
#include "polymake/SparseVector.h"
#include "polymake/SparseMatrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Set.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"

namespace polymake { namespace matroid{

   class Preprocessor{
   private:
      const int interredTrigger;
      const Matrix<Integer> conesFacetValues;
      typedef Map<Vector<Integer>, Set<SparseVector<Integer> > > objmap;
      objmap currentObjects;
      int currentSize;
      int lastInterreduceSize;

      SparseMatrix<Integer> grep_smaller_elements(const Vector<Integer>& facetValues) {
         // Walk through list and collect smaller elements. List is sorted, so
         // we may stop when we reach something bigger.
         return grep_smaller_elements_from_map(facetValues, currentObjects);
      }

      ListMatrix<SparseVector<Integer> > grep_smaller_elements_from_map(const Vector<Integer>& facetValues, const objmap& reductors){
         ListMatrix< SparseVector<Integer> > smaller(0, conesFacetValues.rows());
         Entire<objmap>::const_iterator it = entire(reductors);
         while (!it.at_end() && it->first <= facetValues) {
            cout << "testing " << it->first << " - " << facetValues;
            if (find_if(entire(attach_operation(it->first,facetValues,operations::gt())),operations::non_zero()).at_end()) {
               cout << " smaller!!!111elf" << endl;
               for(Entire< Set< SparseVector<Integer> > >::const_iterator vit = entire(it->second); !vit.at_end(); ++vit){
                  smaller /= *vit;
               }
            } else { cout << endl; }
            ++it;
         }
         return smaller;
      }

      bool is_superfluous(const SparseVector<Integer>& circuit, const Vector<Integer>& facetValues){
         SparseMatrix<Integer> smallerElements = grep_smaller_elements(facetValues);
         int rk = rank(smallerElements);
         int rkExtended = rank(smallerElements/circuit);
         return rk==rkExtended;
      }

      Vector<Integer> intersect_cones(const SparseVector<Integer>& circuit){
         // Intersect the cones belonging to non-zero entries.
         Vector<Integer> maxval(conesFacetValues.cols());
         Matrix<Integer> fvminor(conesFacetValues.minor(indices(circuit),All));
         Entire<Cols<Matrix<Integer> > >::const_iterator colit = entire(cols(fvminor));
         for(Entire<Vector<Integer> >::iterator it=entire(maxval); !it.at_end(); ++it, ++colit) {
            *it = accumulate(*colit, operations::max());
         }
         return maxval;
      }
      


   public:
      Preprocessor(const int& it, const Matrix<Integer>& cfv): interredTrigger(it), conesFacetValues(cfv) {
         lastInterreduceSize = 0;
         currentSize = 0;
      }

      void add(const SparseVector<Integer>& circuit){
         Vector<Integer> facetValues = intersect_cones(circuit);
         bool check = is_superfluous(circuit, facetValues);
         if(!check) {
            currentObjects[facetValues] += circuit;
            currentSize++;
         }
         if(currentSize >= lastInterreduceSize + interredTrigger)
            interreduce();
      }

      void show() {
         cout << currentObjects << endl;
      }

      const objmap& give() const {
         return currentObjects;
      }
      
      void interreduce(){
         // Interreduce all elements and record new size.
         cout << "Interreducing." << endl;
         objmap remember;
         int newSize = 0;
         while(!currentObjects.empty()){
            const Vector<Integer> facetValues(currentObjects.front().first);
            ListMatrix<SparseVector<Integer> > circuits(currentObjects.front().second);
            currentObjects.pop_front();
            ListMatrix<SparseVector<Integer> > smaller = grep_smaller_elements_from_map(facetValues, remember)/grep_smaller_elements_from_map(facetValues, currentObjects);
            Rows<ListMatrix<SparseVector<Integer> > >::iterator cit = rows(circuits).begin();
            while (cit != rows(circuits).end()) {
               SparseVector<Integer> circuit(*cit);
               circuits.delete_row(cit);
               int rk = rank(smaller/circuits);
               int rkExtended = rank(smaller/circuits/circuit);
               if(rk != rkExtended){
                  remember[facetValues] += circuit;
                  newSize++;
                  smaller = smaller/circuit;
               }
               cit = rows(circuits).begin();
            }
         }
         currentObjects = remember;
         lastInterreduceSize = newSize;
         currentSize = newSize;
      }

   };

}
}
#endif
