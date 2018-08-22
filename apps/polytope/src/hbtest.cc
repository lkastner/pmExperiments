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
Matrix<Integer> continued_fraction_to_hilbert_basis(const Vector<Integer>& dcf);
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
Array<Vector<Rational> > cone_to_divisor_slices(const Matrix<Rational>& rays, const Matrix<Rational>& facets);
Vector<Rational> get_oneStep_slice(const Matrix<Rational>& rays, const Matrix<Rational>& facets);
Integer find_number_of_steps(const Matrix<Rational>& rays, const Matrix<Rational>& facets);
Vector<Integer> continued_fraction_from_rational(const Rational& r);
std::pair<Integer, Vector<Integer> > find_index_of_divisor(const Integer& ray0val, const Integer& ray1val, const Integer& n, const Integer& q);
Matrix<Integer> project_and_lift_3dim(const Array<Vector<Rational> >& slices, const Matrix<Integer>& rays, const Matrix<Integer>& facets, int facetIndex);
Matrix<Integer> threedim_hb_old(const Matrix<Integer>& rays, const Matrix<Integer>& facets);
Matrix<Integer> threedim_hb_v1(const Matrix<Integer>& rays, const Matrix<Integer>& facets);
bool check_containment(const Matrix<Integer>& reductor, const Vector<Integer>& a, const Matrix<Integer>& facets);
Matrix<Integer> reduce(const Matrix<Integer>& reductor, const Matrix<Integer>& given, const Matrix<Integer>& facets);
int resp_containment(const Vector<Integer>& a, const Vector<Integer>& b, const Matrix<Integer>& facets);
Matrix<Integer> interred(const Matrix<Integer>& gens, const Matrix<Integer>& facets);
Integer optimize_upper_slice_bound(const Matrix<Integer>& gens, const Matrix<Integer>& facets);


/////////////////////////////////////////////////
// Global variables
//
Map<Vector<Integer>, Matrix<Integer> > cqshash;


/////////////////////////////////////////////////
// Classes
//
class twofacet{
private:
   Vector<Integer> facet, dcf;
   std::pair<SparseMatrix<Integer>, SparseMatrix<Integer> > orthBasis;
   std::pair<Matrix<Integer>, Matrix<Integer> > twodimstd;
   bool smooth;
   Integer n, q;
   Matrix<Integer> projection, generatingRays, twodimFacets, hb;

public:
   twofacet(const Matrix<Integer>& rays, const Matrix<Integer>& facets, int facetIndex){
      int rayIndex;
      facet = facets.row(facetIndex);
      orthBasis = find_orthogonal_lattice_basis(facet);
      for(rayIndex = 0; rayIndex < 3; rayIndex++){
         if(facet * rays.row(rayIndex) != 0){
            break;
         }
      }
      generatingRays = rays.minor(~scalar2set(rayIndex), All);
      // cout << "GR: " << endl << generatingRays << endl;
      generatingRays = T(orthBasis.first).minor(~scalar2set(0),All) * T(generatingRays);
      // cout << "GR: " << endl << generatingRays << endl;
      smooth = (abs(det(generatingRays)) == 1);
      twodimstd = twodim_standard_form(generatingRays.col(0), generatingRays.col(1));
      projection = twodimstd.second * T(orthBasis.first).minor(~scalar2set(0),All);
      generatingRays = twodimstd.second * generatingRays;
      // cout << "GR: " << endl << generatingRays << endl;
      n = generatingRays(0,1);
      q = generatingRays(1,1);
      twodimFacets = Matrix<Integer>(2,2);
      twodimFacets(0,0) = 1;
      twodimFacets(0,1) = 0;
      twodimFacets(1,0) = -q;
      twodimFacets(1,1) = n;
      if(!smooth){
         dcf = continued_fraction_from_rational(Rational(n,n-q));
         hb = continued_fraction_to_hilbert_basis(dcf);
      } else {
         hb = generatingRays;
      }
      // cout << "n: " << n << " q: " << q << endl;
   }

   Matrix<Integer> get_lift_of_pts(const Matrix<Integer>& dim2pts, Integer ht){
      Matrix<Integer> result(dim2pts);
      result = result * T(twodimstd.first);
      result = (ht * ones_vector<Integer>(result.rows())) | result;
      result = result * orthBasis.second;
      return result;
   }

   Matrix<Integer> get_facet_hb(){
      return get_lift_of_pts(hb, 0);
   }

   Matrix<Integer> get_gens_of_slice(const Vector<Rational>& slice){
      // Vars:
      Vector<Rational> vertex, divisor;
      Vector<Integer> facetVals;
      Matrix<Integer> twodimGens;
      Rational height;
      std::pair<Integer, Vector<Integer> > twodimIndex;
      // Proc:
      height = facet * slice;
      vertex = projection * slice;
      divisor = find_biggest_contained_divisor(T(generatingRays), twodimFacets, vertex);
      if(smooth or ((denominator(divisor[0]) == 1) and ((denominator(divisor[1]) == 1)))){
         twodimGens = Matrix<Integer>(1,2);
         twodimGens(0,0) = numerator(divisor[0]);
         twodimGens(0,1) = numerator(divisor[1]);
      } else {
         // cout << vertex << " - " << divisor << endl;
         // cout << twodimFacets*vertex << " - " << twodimFacets*divisor << endl;
         facetVals = twodimFacets*divisor;
         twodimIndex = find_index_of_divisor(facetVals[0], facetVals[1], n, q);
         // cout << "has index: " << twodimIndex.first << endl;
         twodimGens = compute_generators_of_ei(twodimIndex.first, dcf, hb);
         twodimGens -= repeat_row(twodimIndex.second, twodimGens.rows());
      }
      // cout << "Gens new: " << endl << twodimGens << endl;
      // cout << "Eval: " << endl << twodimFacets*T(twodimGens) << endl;
      return get_lift_of_pts(twodimGens, numerator(height));
   }

};


/////////////////////////////////////////////////
// Implementation
//
Matrix<Integer> threedim_hb_v1(const Matrix<Integer>& rays, const Matrix<Integer>& facets){
   // cout << "V1 entered." << endl;
   // cout << "Rays: " << endl << rays << endl;
   // cout << "Facets: " << endl << facets << endl;
   Matrix<Rational> raysR(rays), facetsR(facets);
   Matrix<Integer> result, newGens;
   Vector<Rational> oneStep = get_oneStep_slice(raysR, facetsR);
   Integer bound = find_number_of_steps(raysR, facetsR), current = 1, tmp;
   twofacet t0(rays, facets, 0), t1(rays, facets, 1), t2(rays, facets, 2);
   // cout << "Init done." << endl;
   // Initializing w/ facethbs:
   result = t0.get_facet_hb();
   // cout << "First hb: " << endl << t0.get_facet_hb() << endl;
   result = result / t1.get_facet_hb();
   // cout << "Second hb: " << endl << t1.get_facet_hb() << endl;
   result = result / t2.get_facet_hb();
   // cout << "Third hb: " << endl << t2.get_facet_hb() << endl;
   // cout << "Init HB size: " << result.rows() << endl;
   // cout << result << endl;
   result = interred(result, facets);
   // cout << "Reduced HB size: " << result.rows() << endl;
   // cout << result << endl;
   bound = optimize_upper_slice_bound(result, facets);
   // cout << "Setup done." << endl;
   // Vile loop:
   while(current <= bound){
      newGens = t0.get_gens_of_slice(current * oneStep);
      newGens = newGens / t1.get_gens_of_slice(current * oneStep);
      newGens = newGens / t2.get_gens_of_slice(current * oneStep);
      // cout << "Computed gens." << endl;
      newGens = reduce(result, newGens, facets);
      newGens = interred(newGens, facets);
      if(newGens.rows() > 0){
         tmp = optimize_upper_slice_bound(newGens, facets);
         // cout << "Bound optimized." << endl;
         bound = bound < tmp ? bound : tmp;
      }
      result = result / newGens;
      current++;
      cout << "Current: " << current << " bound: " << bound << endl;
   }
   // cout << "Vile loop done." << endl;
   return result;
}


Integer optimize_upper_slice_bound(const Matrix<Integer>& gens, const Matrix<Integer>& facets){
   Integer max, result(accumulate(facets * gens.row(0), operations::max()));
   for(const auto& g : rows(gens)) {
      max = accumulate(facets * g, operations::max());
      result = result < max ? result : max;
   }
   return result;
}


Matrix<Integer> interred(const Matrix<Integer>& gens, const Matrix<Integer>& facets){
   Matrix<Integer> result(0, gens.rows()), keep(0, gens.rows()), reductor(gens);
   Vector<Integer> current;
   int check;
   bool bad;
   while(0 < reductor.rows()){
      current = reductor.row(0);
      reductor = reductor.minor(~scalar2set(0), All);
      keep = Matrix<Integer>(0, gens.rows());
      bad = false;
      for(const auto& r : rows(reductor)) {
         check = resp_containment(current, r, facets);
         if(check == 1){
            bad = true;
            keep = keep / r;
         } else if(check == 0){
            keep = keep / r;
         }
      }
      if(!bad) result = result / current;
      reductor = keep;
      // cout << reductor.rows() << endl;
   }
   return result;
}


int resp_containment(const Vector<Integer>& a, const Vector<Integer>& b, const Matrix<Integer>& facets){
   // Output:
   // 1   if  a\in b+\sigma
   // -1  if  b\in a+\sigma
   // 0   else
   Vector<Integer> test = facets*(a-b);
   Integer min, max;
   min = accumulate(test, operations::min());
   max = accumulate(test, operations::max());
   if(min >= 0) return 1;
   if(max <= 0) return -1;
   return 0;
}


Matrix<Integer> threedim_hb_old(const Matrix<Integer>& rays, const Matrix<Integer>& facets){
   Matrix<Rational> raysR(rays), facetsR(facets);
   Matrix<Integer> result, current;
   cout << "Slicing." << endl;
   Array<Vector<Rational> > slices(cone_to_divisor_slices(raysR, facetsR));
   cout << slices.size() << endl;
   result = project_and_lift_3dim(slices, rays, facets, 0);
   cout << "Ok." << endl;
   cout << "Upper slice bound could be: " << optimize_upper_slice_bound(result, facets) << endl;
   current = project_and_lift_3dim(slices, rays, facets, 1);
   //current = reduce(result, current, facetsR);
   result = result / current;
   cout << "Ok." << endl;
   cout << "Upper slice bound could be: " << optimize_upper_slice_bound(result, facets) << endl;
   current =  project_and_lift_3dim(slices, rays, facets, 2);
   //current = reduce(result, current, facetsR);
   result = result / current;
   cout << "Ok." << endl;
   return interred(result, facets);
}



Matrix<Integer> project_and_lift_3dim(const Array<Vector<Rational> >& slices, const Matrix<Integer>& rays, const Matrix<Integer>& facets, int facetIndex){
   cout << "Entered pal3dim." << endl;
   twofacet testing(rays, facets, facetIndex);
   Matrix<Integer> twodimGens, reductor;
   Vector<Integer> facet(facets.row(facetIndex));
   Rational height;

   reductor = testing.get_facet_hb();
   cout << "Setup done." << endl;
   for(const auto& slice : slices) {
      height = facet * slice;
      // cout << *slice << " at height " << height << endl;
      if(height == 0){
         continue;
      } else {
         twodimGens = testing.get_gens_of_slice(slice);
         twodimGens = reduce(reductor, twodimGens, facets);
         reductor = reductor / twodimGens;
      }
   }
   return reductor;
}


Matrix<Integer> reduce(const Matrix<Integer>& reductor, const Matrix<Integer>& given, const Matrix<Integer>& facets){
   bool test;
   Matrix<Integer> result(0,reductor.cols());
   for(const auto& gen : rows(given)) {
      test = check_containment(reductor, gen, facets);
      if(!test){
         result = result/(Vector<Integer>(gen));
      }
   }
   return result;
}


bool check_containment(const Matrix<Integer>& reductor, const Vector<Integer>& a, const Matrix<Integer>& facets){
   Vector<Rational> test;
   Rational min;
   for(const auto& r : rows(reductor)) {
      test = facets*(a-(r));
      // Should just be min(test)
      min = accumulate(test, operations::min());
      if(min >= 0) return true;
   }
   return false;
}


std::pair<Integer, Vector<Integer> > find_index_of_divisor(const Integer& ray0val, const Integer& ray1val, const Integer& n, const Integer& q){
   ExtGCD<Integer> gcd = ext_gcd(q, n);
   Integer result = ray0val + gcd.p*ray1val;
   Vector<Integer> shift(2);
   shift[0] = gcd.p*ray1val;
   shift[1] = (shift[0]*q-ray1val)/n;
   return std::pair<Integer, Vector<Integer> >(result, shift);
}


Vector<Integer> continued_fraction_from_rational(const Rational& r){
   Integer n(numerator(r)), q(denominator(r)), c(ceil(r));
   Rational tmp;
   if(q==1){
      Vector<Integer> result(1);
      result[0] = n;
      return result;
   } else {
      tmp = c - r;
      return c | continued_fraction_from_rational(1/tmp);
   }
}


Integer find_number_of_steps(const Matrix<Rational>& rays, const Matrix<Rational>& facets){
   Matrix<Rational> prod = facets * (T(rays));
   Rational result(prod(0,0)), temp;
   int i, j;
   for(i=0; i<prod.rows(); i++){
      result = result > prod(i,0) ? result : prod(i,0);
   }
   for(i=0; i<prod.rows(); i++){
      temp = 0;
      for(j=0; j<prod.cols(); j++){
         temp = temp > prod(i,j) ? temp : prod(i,j);
      }
      result = result < temp ? result : temp;
   }
   return floor(result);
}


Array<Vector<Rational> > cone_to_divisor_slices(const Matrix<Rational>& rays, const Matrix<Rational>& facets){
   int bound = (int) (find_number_of_steps(rays, facets));
   Array<Vector<Rational> > result(bound);
   Vector<Rational> oneStep = get_oneStep_slice(rays, facets);
   for(int i = 0; i < bound; i++){
      result[i] = i*oneStep;
   }
   return result;
}


Vector<Rational> get_oneStep_slice(const Matrix<Rational>& rays, const Matrix<Rational>& facets){
   Rational test;
   int n = rays.cols();
   Vector<Rational> oneStep(n);
   for(const auto& facet : rows(facets)) {
      for(const auto& ray : rows(rays)) {
         test = ray * facet;
         if(test > 0){
            oneStep += (1/test) * ray;
         }
      }
   }
   return oneStep;
}


Vector<Rational> find_biggest_contained_divisor(const Matrix<Integer>& rays, const Matrix<Integer>& facets, const Vector<Rational> input){
   Vector<Rational> result(input), facetVals;
   Rational test, fix;
   Map<Vector<Integer>, Vector<Integer> > facetRay;
   for(const auto& facet : rows(facets)) {
      for(const auto& ray : rows(rays)) {
         test = ray * facet;
         if(test > 0){
            facetRay[facet] = ray;
         } 
      }
   }
   facetVals = facets * input;
   for(const auto& facet : rows(facets)) {
      test = input * facet;
      if(denominator(test) != 1){
         //cout << "This is not an int." << endl;
         fix = (ceil(test) - test) / (facetRay[facet] * facet);
         result += fix * facetRay[facet];
      }
   }
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
   Matrix<Integer> MR(unit_matrix<Integer>(2)), ML(unit_matrix<Integer>(2)),  changeCoord(2,2), mirror(2,2), add(2,2), sub(2,2);
   make_matrix_pair(L, R, gen0, 1);
   Integer mod, result, factor;
   MR.multiply_from_right(R);
   ML.multiply_from_right(L);
   Integer n;
   Vector<Integer> imGen0(gen0), imGen1(gen1), test;
   changeCoord(0,0) = 0;
   changeCoord(0,1) = 1;
   changeCoord(1,0) = 1;
   changeCoord(1,1) = 0;
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
   result = imGen1[1];
   mod = result % n;
   if(mod < 0) mod += n;
   factor = (mod - result) / n;
   add(0,0) = 1;
   add(0,1) = 0;
   add(1,0) = factor;
   add(1,1) = 1;
   sub(0,0) = 1;
   sub(0,1) = 0;
   sub(1,0) = -factor;
   sub(1,1) = 1;
   test = add * imGen1;
   MR = add * MR;
   ML = ML * sub;
   // cout << "Testing: " << imGen1 << " -- " << test << endl;
   return std::pair<Matrix<Integer>, Matrix<Integer> >(ML,MR);
}


Matrix<Integer> compute_generators_of_ei(Integer i, const Vector<Integer>& dcf, const Matrix<Integer>& hb){
   // cout << "canonical: " << canonical << endl;
   Vector<Integer> last(hb.row(hb.rows() - 1)), key;
   Integer n(last[0]);
   Matrix<Integer> result;
   std::pair<Integer, Vector<Integer> > fix(adjust_entry(i, last));
   if((fix.first == n) or (fix.first == 0)) cout << "Take care!" << endl;
   key = Vector<Integer>(fix.first | dcf);
   if(cqshash.exists(key)){
      result = cqshash[key];
   } else {
      result = compute_ext_degrees(fix.first, -1, dcf, hb) / hb.row(hb.rows() - 1);
      cqshash[key] = result;
   }
   return result + repeat_row(fix.second, result.rows());
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
   Integer result = i, n(last[0]), mod, factor;
   mod = result % n;
   if(mod <= 0) mod += n;
   factor = (mod - result) / n;
   result += factor * n;
   shift -= factor * last;
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
      newdcf = dcf.slice(range(0,size-1));
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

Function4perl(&cone_to_divisor_slices, "cone_to_divisor_slices");

Function4perl(&continued_fraction_from_rational, "continued_fraction_from_rational");

Function4perl(&threedim_hb_old, "threedim_hb_old");

Function4perl(&threedim_hb_v1, "threedim_hb_v1");

} // namespace polymake
} // namespace polytope
