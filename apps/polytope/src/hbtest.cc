/* Copyright (c) 1997-2015
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

#include "polymake/client.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/vector"
#include "polymake/PowerSet.h"
#include "polymake/common/lattice_tools.h"


namespace polymake {
namespace polytope {
namespace {

Matrix<Integer> module_generator_recursion(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
Matrix<Integer> eq2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
Matrix<Integer> gt2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);

Matrix<Integer> module_generator_recursion(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   cout << "Entering recursion." << dcf << " i: " << i << " j: " << j << endl;
   int size = dcf.size();
   if(size==0){
      cout << "Returning empty mat." << endl;
      return Matrix<Integer>(0,2);
   }
   if(dcf[size-1] == 2) return eq2_modgen(i, j, dcf, hb);
   else if(dcf[size-1] > 2) return gt2_modgen(i, j, dcf, hb);
   else throw 42;
}

Matrix<Integer> eq2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   int size = dcf.size();
   cout << "==2 case." << endl;
   Integer n, q, ntilda, qtilda, diff, newi;
   Vector<Integer> last(hb.row(size+1)), newlast(hb.row(size)), newdcf, fix;
   Matrix<Integer> newhb, result, fixmat;
   n = last[0];
   q = last[1];
   ntilda = newlast[0];
   qtilda = newlast[1];
   if((i <= ntilda) && (j <= ntilda)){
      newdcf = dcf.slice(0,size-1);
      if(size == 1){
         newdcf = Vector<Integer>(0);
      }
      cout << "dcf: " << dcf << " newdcf: " << newdcf << endl;
      newhb = hb.minor(~scalar2set(size+1), All);
      result = module_generator_recursion(i, j, newdcf, newhb);
      if( i+j <= n) result /= last;
      return result;
   } else if (ntilda < i){
      diff = n - ntilda;
      fix = last - newlast;
      newi = i - diff;
      result = module_generator_recursion(newi, j, dcf, hb);
      if(i+j > n){
         cout << "Previous result: " << endl << result << endl << "------------" << endl;
         cout << "Last? row: " << result.row(result.rows() - 1) << endl;
         cout << "n and q: " << n << " " << q << endl;
         cout << "i and j: " << i << j << " dcf: " << dcf << endl;
         result = result.minor(~scalar2set(result.rows() - 1), All);
      }
      fixmat = repeat_row(fix, result.rows());
      result += fixmat;
      return result;
   } else {
      return module_generator_recursion(j, i, dcf, hb);
   }
}

Matrix<Integer> gt2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   int size = dcf.size();
   Integer n, q, ntilda, qtilda, diff, newi, newj;
   Vector<Integer> last(hb.row(size+1)), secondlast(hb.row(size)), newlast, newdcf;
   Matrix<Integer> result, newhb;
   newlast = last - secondlast;
   n = last[0];
   q = last[1];
   ntilda = newlast[0];
   qtilda = newlast[1];
   diff = n - ntilda;
   if((diff < i) && (diff < j)){
      newi = i - diff;
      newj = j - diff;
      newdcf = Vector<Integer>(dcf);
      newdcf[size-1]--;
      newhb = hb.minor(~scalar2set(size+1), All) / newlast;
      result = module_generator_recursion(newi, newj, newdcf, newhb);
      if(result.rows() > 0){
         result += repeat_row(2*secondlast, result.rows());
      }
      if(i+j <= n){
         result /= last;
      }
      return result;
   } else if (i <= diff) {
      newi = i + ntilda;
      result = module_generator_recursion(newi, j, dcf, hb);
      result += repeat_row(secondlast, result.rows());
      if(i+j < n){
         result /= last;
      }
      return result;
   } else {
      return module_generator_recursion(j, i, dcf, hb);
   }
}

Matrix<Integer> continued_fraction_to_hilbert_basis(const Vector<Integer>& dcf){
   int n = dcf.size(), i;
   Matrix<Integer> hb(n+2, 2);
   hb(0,1) = 1;
   hb(1,0) = 1;
   hb(1,1) = 1;
   for(i=0; i<n; i++){
      hb(i+2,0) = dcf[i] * hb(i+1,0) - hb(i,0);
      hb(i+2,1) = dcf[i] * hb(i+1,1) - hb(i,1);
   }
   return hb;
}

} // namespace

Function4perl(&continued_fraction_to_hilbert_basis, "continued_fraction_to_hilbert_basis");

Function4perl(&module_generator_recursion, "module_generator_recursion");

} // namespace polymake
} // namespace polytope
