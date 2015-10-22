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

class HBStep {
public:
   Array<std::pair<Vector<Rational>, Vector<Rational> > > modules;
   Matrix<Integer> rays;
   Matrix<Integer> facets;


   
};

/////////////////////////////////////////////////
// Declaration
//
Matrix<Integer> module_generator_recursion(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
Matrix<Integer> eq2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
Matrix<Integer> gt2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
std::pair<Integer, Vector<Integer> > adjust_entry(Integer i, const Vector<Integer>& last);
std::pair<std::pair<Integer, Integer>, Vector<Integer> > adjust_entries(Integer i, Integer j, const Matrix<Integer>& hb);
std::pair<Integer, Vector<Integer> > get_canonical_divisor(const Matrix<Integer>& hb);
Matrix<Integer> compute_ext_degrees(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb); 
Matrix<Integer> compute_generators_of_ei(Integer i, const Vector<Integer>& dcf, const Matrix<Integer>& hb);
std::pair<Matrix<Integer>, Matrix<Integer> > twodim_standard_form(const Vector<Integer>& gen1, const Vector<Integer>& gen2);
Integer make_matrix_pair(SparseMatrix2x2<Integer>& R, SparseMatrix2x2<Integer>& L, const Vector<Integer>& u, int i);
std::pair<SparseMatrix<Integer>, SparseMatrix<Integer> > find_orthogonal_lattice_basis(const Vector<Integer>& primitive);
Vector<Rational> find_biggest_contained_divisor(const Matrix<Integer>& rays, const Matrix<Integer>& facets, const Vector<Rational> input);

/////////////////////////////////////////////////
// Implementation
//
Vector<Rational> find_biggest_contained_divisor(const Matrix<Integer>& rays, const Matrix<Integer>& facets, const Vector<Rational> input){
   Vector<Rational> result(input);

   return result;
}

std::pair<SparseMatrix<Integer>, SparseMatrix<Integer> > find_orthogonal_lattice_basis(const Vector<Integer>& primitive){
   Integer g;
   Vector<Integer> currentVec(primitive);
   int i, dim = primitive.size();
   SparseMatrix<Integer> right = unit_matrix<Integer>(dim), left = unit_matrix<Integer>(dim);
   SparseMatrix2x2<Integer> R, L;
   for(i=1; i<dim; i++){
      if(primitive[i] == 0){
         continue;
      }
      g = make_matrix_pair(L, R, currentVec, i);
      right.multiply_from_left(R);
      left.multiply_from_right(L);
      currentVec[0] = g;
   }
   return std::pair<SparseMatrix<Integer>, SparseMatrix<Integer> >(left, right);
}



Integer make_matrix_pair(SparseMatrix2x2<Integer>& L, SparseMatrix2x2<Integer>& R, const Vector<Integer>& u, int i){
   ExtGCD<Integer> gcd = ext_gcd(u[0], u[i]);
   R.i = 0;
   R.j = i;
   R.a_ii = gcd.p;
   R.a_ij = gcd.q;
   R.a_ji = -gcd.k2;
   R.a_jj = gcd.k1;
   L.i = 0;
   L.j = i;
   L.a_ii = gcd.k1;
   L.a_ij = -gcd.q;
   L.a_ji = gcd.k2;
   L.a_jj = gcd.p;
   return gcd.g;
}


std::pair<Matrix<Integer>, Matrix<Integer> > twodim_standard_form(const Vector<Integer>& gen0, const Vector<Integer>& gen1){
   SparseMatrix2x2<Integer> R, L;
   Matrix<Integer> MR(unit_matrix<Integer>(2)), ML(unit_matrix<Integer>(2)), addone(2,2), subone(2,2), changeCoord(2,2), mirror(2,2);
   make_matrix_pair(L, R, gen0, 1);
   MR.multiply_from_right(R);
   ML.multiply_from_right(L);
   Integer n;
   Vector<Integer> imGen0(gen0), imGen1(gen1);
   changeCoord(0,0) = 0;
   changeCoord(0,1) = 1;
   changeCoord(1,0) = 1;
   changeCoord(1,1) = 0;
   addone(0,0) = 1;
   addone(0,1) = 0;
   addone(1,0) = 1;
   addone(1,1) = 1;
   subone(0,0) = 1;
   subone(0,1) = 0;
   subone(1,0) = -1;
   subone(1,1) = 1;
   ML = ML * changeCoord;
   MR = changeCoord * MR;
   imGen1 = MR * imGen1;
   if(imGen1[0] < 0){
      mirror(0,0) = -1;
      mirror(0,1) = 0;
      mirror(1,0) = 0;
      mirror(1,1) = 1;
      ML = ML * mirror;
      MR = mirror * MR;
      imGen1 = mirror * imGen1;
   }
   n = imGen1[0];
   while(imGen1[1] < 0){
      ML = ML * subone;
      MR = addone * MR;
      imGen1 = addone * imGen1;
   }
   while(imGen1[1] >= n){
      ML = ML * addone;
      MR = subone * MR;
      imGen1 = subone * imGen1;
   }
   return std::pair<Matrix<Integer>, Matrix<Integer> >(ML,MR);
}


Matrix<Integer> compute_generators_of_ei(Integer i, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   // cout << "canonical: " << canonical << endl;
   return compute_ext_degrees(i, -1, dcf, hb) / hb.row(hb.rows() - 1);
}


Matrix<Integer> compute_ext_degrees(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   Integer newi, cj, newj;
   Matrix<Integer> result;
   std::pair<Integer, Vector<Integer> > canonical = get_canonical_divisor(hb);
   cj = canonical.first + j;
   std::pair<std::pair<Integer, Integer>, Vector<Integer> > adjustment = adjust_entries(i, cj, hb);
   // cout << "Adjusted: " << adjustment << endl;
   newi = adjustment.first.first;
   newj = adjustment.first.second;
   result = module_generator_recursion(newi, newj, dcf, hb);
   result += repeat_row(canonical.second + adjustment.second, result.rows());
   return result;
}



std::pair<Integer, Vector<Integer> > get_canonical_divisor(const Matrix<Integer>& hb){
   Vector<Integer> last(hb.row(hb.rows() -1)), secondlast(hb.row(hb.rows() -2)), shift(2);
   Integer n(last[0]), q(last[1]), index;
   index = n - secondlast[0] + 1;
   shift[0] = 1-index;
   shift[1] = (q+1-index*q)/n;
   return std::pair<Integer, Vector<Integer> >(index, shift);
}

std::pair<std::pair<Integer, Integer>, Vector<Integer> > adjust_entries(Integer i, Integer j, const Matrix<Integer>& hb){
   Vector<Integer> last(hb.row(hb.rows()-1)), shift(2);
   Integer n(last[0]), newi, newj;
   std::pair<Integer, Vector<Integer> > adjust = adjust_entry(i, last);
   // cout << "Adjusted: " << adjust << endl;
   newi = adjust.first;
   shift += adjust.second;
   adjust = adjust_entry(j, last);
   // cout << "Adjusted: " << adjust << endl;
   newj = adjust.first;
   shift += adjust.second;
   return std::pair<std::pair<Integer, Integer>, Vector<Integer> >(std::pair<Integer, Integer>(newi, newj), shift);
}


std::pair<Integer, Vector<Integer> > adjust_entry(Integer i, const Vector<Integer>& last){
   Vector<Integer> shift(2);
   Integer result = i, n(last[0]);
   while(result < 1){
      result += n;
      shift -= last;
   }
   while(result >= n){
      result -= n;
      shift += last;
   }
   return std::pair<Integer, Vector<Integer> >(result, shift);
}


Matrix<Integer> module_generator_recursion(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   // cout << "Entering recursion." << dcf << " i: " << i << " j: " << j << endl;
   int size = dcf.size();
   if(size==0){
      // cout << "Returning empty mat." << endl;
      return Matrix<Integer>(0,2);
   }
   if(dcf[size-1] == 2) return eq2_modgen(i, j, dcf, hb);
   else if(dcf[size-1] > 2) return gt2_modgen(i, j, dcf, hb);
   else throw 42;
}

Matrix<Integer> eq2_modgen(Integer i, Integer j, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   int size = dcf.size();
   // cout << "==2 case. dcf: " << dcf << endl;
   Integer n, q, ntilda, qtilda, diff, newi;
   Vector<Integer> last(hb.row(size+1)), newlast(hb.row(size)), newdcf, fix;
   // cout << "last: " << last << " newlast: " << newlast << endl;
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
      // cout << "dcf: " << dcf << " newdcf: " << newdcf << endl;
      newhb = hb.minor(~scalar2set(size+1), All);
      result = module_generator_recursion(i, j, newdcf, newhb);
      if( i+j <= n) result /= last;
      return result;
   } else if (ntilda < i){
      diff = n - ntilda;
      fix = last - newlast;
      newi = i - diff;
      result = module_generator_recursion(newi, j, dcf, hb);
      if(newi+j <= n){
         // cout << "Previous result: " << endl << result << endl << "------------" << endl;
         // cout << "Last? row: " << result.row(result.rows() - 1) << endl;
         // cout << "n and q: " << n << " " << q << endl;
         // cout << "i and j: " << i << j << " dcf: " << dcf << endl;
         if((result.rows() == 0) || (n != (result.row(result.rows() - 1))[0])){
            cout << "Something went wrong!" << endl;
         }
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
      result -= repeat_row(newlast, result.rows());
      if(i+j <= n){
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

Function4perl(&compute_ext_degrees, "compute_ext_degrees");

Function4perl(&compute_generators_of_ei, "compute_generators_of_ei");

Function4perl(&twodim_standard_form, "twodim_standard_form");

} // namespace polymake
} // namespace polytope
