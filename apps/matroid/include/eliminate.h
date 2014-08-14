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
      const Matrix<Integer> facets, conesFacetValues;
      typedef Map<Vector<Integer>, Set<SparseVector<Integer> > > objmap;
      objmap currentObjects;
      int lastInterreduceSize;

      SparseMatrix<Integer> grep_smaller_elements(const SparseVector<Integer>& circuit, const Vector<Integer>& facetValues) {
         // Walk through list and collect smaller elements. List is sorted, so
         // we may stop when we reach something bigger.
         ListMatrix< SparseVector<Integer> > smaller(0, conesFacetValues.rows());
         Entire<objmap>::const_iterator it = entire(currentObjects);
         while (!it.at_end() && it->first > facetValues) {
            if (/*componentval less first - facetvalues*/1) {
               for(Entire< Set< SparseVector<Integer> > >::const_iterator vit = entire(it->second); !vit.at_end(); ++vit){
                  smaller /= *vit;
               }
            }
            ++it;
         }
         return smaller;
      }

      bool is_superflous(const SparseVector<Integer>& circuit, const Vector<Integer>& facetValues){
         SparseMatrix<Integer> smallerElements = grep_smaller_elements(circuit, facetValues);
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
      
      void interreduce(){
         // Interreduce all elements and record new size.

      }


   public:
      Preprocessor(const int& it, const Matrix<Integer>& f, const Matrix<Integer>& cfv): interredTrigger(it), facets(f), conesFacetValues(cfv) { }

      void add(const SparseVector<Integer>& circuit){
         Vector<Integer> facetValues = intersect_cones(circuit);
         bool check = is_superflous(circuit, facetValues);
         if(!check) 
            currentObjects[facetValues] += circuit;
         if(currentObjects.size() >= lastInterreduceSize + interredTrigger)
            interreduce();
      }

      void show() {
         cout << currentObjects << endl;
      }

   };

}
}
#endif
