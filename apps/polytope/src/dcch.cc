/* Copyright (c) 2021
   Lars Kastner

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
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Map.h"
#include "polymake/linalg.h"
#include "polymake/common/lattice_tools.h"
#include "polymake/polytope/convex_hull.h"

#include <random>
#include <queue>



namespace polymake {
namespace polytope {
namespace {

template<typename Scalar>
void normalize_facet_vector(Vector<Scalar>& facet){
   Int i = 0;
   while(i<facet.dim() && facet[i] == 0) i++;
   if(i == facet.dim()) throw std::runtime_error("Cannot normalize zero vector!");
   if(facet[i] > 0) facet *= 1/facet[i];
   else if(facet[i] < 0) facet *= -1/facet[i];
}

using polymake::common::primitive;

class DCCH_Logger {
   private:
      Int depth;
      Int level;

   public:
      DCCH_Logger() : depth(0), level(0) {}
      DCCH_Logger(Int a) : depth(0), level(a) {}

      void descend() {
         depth++;
      }

      void backtrack() {
         depth--;
      }

      void print_indent() {
         for(int i=0; i<depth; i++){
            cout << "| ";
         }
      }

      void log(std::string text, Int l) {
         if(l>=level) return;
         print_indent();
         cout << text << endl;
      }

      template<typename Scalar>
      void log(const Vector<Scalar>& v, Int l){
         if(l>=level) return;
         print_indent();
         for(const auto& entry : v){
            cout << entry << " ";
         }
         cout << endl;
      }

};

template<typename Scalar>
void orthogonalize(Vector<Scalar>& result, const Matrix<Scalar>& orth_affine_hull){
   for(const auto& row : rows(orth_affine_hull)){
      Scalar rowval = row*row;
      if(rowval != 0){
         result -= ((result*row)/rowval) * row;
      }
   }
}


template<typename Scalar>
class LiftingStack {
   private:
      using keytype = std::tuple<const Matrix<Scalar>&, const Matrix<Scalar>&, const Vector<Scalar>&>;
      std::vector<keytype> stack;

   public:
      LiftingStack() {}

      void push(const Matrix<Scalar>& pts, const Matrix<Scalar>& oah, const Vector<Scalar>& facet){
         keytype t(pts, oah, facet);
         stack.push_back(t);
      }

      void pop() {
         stack.pop_back();
      }


      Vector<Scalar> lift_facet(const Vector<Scalar>& facetFacet){
         Vector<Scalar> result(facetFacet);
         for(Int i=0; i<(int)stack.size(); i++){
            lift_facet(result, stack.size()-1-i);
         }
         return result;
      }
      
      void lift_facet(Vector<Scalar>& result, Int i){
         const auto& [points, orth_affine_hull, facet] = stack[i];
         Scalar factor(0);
         bool found = false;
         for(const auto& pt : rows(points)){
            Scalar fval = pt*facet;
            if(fval != 0){
               Scalar check = -(result*pt) / fval;
               if(!found){
                  factor = check;
                  found = true;
               }
               else if(check > factor){
                  factor = check;
               }
            }
         }
         result += factor*facet;
         orthogonalize(result, orth_affine_hull);
         normalize_facet_vector(result);
      }

};

template<typename Scalar>
struct MinMax {
   public:
      Scalar min, max;
      Int min_pt, max_pt;
};

template<typename Scalar>
class DCCH {
   private:
      const Matrix<Scalar>& points;
      const Vector<Scalar>& positive_hyperplane;
      const Int threshold;
      DCCH_Logger& logger;
      LiftingStack<Scalar>& ls;
      Int dim;
      Set<Int> affine_hull;
      Matrix<Scalar> orth_affine_hull;
      Set<Vector<Scalar>> facets, known;
      std::queue<Vector<Scalar>> queue;


      
      MinMax<Scalar> find_min_max(const Vector<Scalar>& h){
         logger.log("find_min_max", 3);
         MinMax<Scalar> result;
         result.max = h * points.row(0);
         result.min = h * points.row(0);
         result.max_pt = 0;
         result.min_pt = 0;
         for(Int i=0; i<points.rows(); i++){
            Scalar check = h * points.row(i);
            if(check > result.max){
               result.max = check;
               result.max_pt = i;
            }
            if(check < result.min){
               result.min = check;
               result.min_pt = i;
            }
         }
         logger.log("find_min_max done", 3);
         return result;
      }


      Scalar find_tilting_factor(const Vector<Scalar>& hyp, const Vector<Scalar>& tilt){
         logger.log("find_tilting_factor", 3);
         // TODO: Why?
         MinMax<Scalar> tiltmm(find_min_max(tilt));
         Vector<Scalar> tilttmp(tilt);
         if(tiltmm.max < 0){
            tilttmp *= -1;
         }
         bool result_set = false;
         Scalar result(0);
         for(const auto& pt : rows(points)){
            Scalar denom = pt * tilttmp;
            // Only use positive side.
            if(denom > 0){
               Scalar factor = - (pt*hyp) / denom;
               if(factor > result) result = factor;
               if(!result_set){
                  result = factor;
                  result_set = true;
               }
            }
         }

         if(tiltmm.max < 0){
            result *= -1;
         }
         logger.log("find_tilting_factor done", 3);
         return result;
      }


      Vector<Scalar> find_next_basis_vector(const Vector<Scalar>& prev, const ListMatrix<Vector<Scalar>>& facet_affine_hull){
         Matrix<Scalar> gens = null_space(prev / facet_affine_hull);
         Vector<Scalar> result(zero_vector<Scalar>(gens.cols()));
         if(gens.rows() == 0){
            throw std::runtime_error("No next basis vector can be found!");
         }
         while(is_zero(result)){
            for(const auto& gen : rows(gens)){
               // TODO: Magic numbers
               result += ((rand() % 10) - 5) * gen;
            }
         }
         MinMax<Scalar> resultmm(find_min_max(result));
         return result;
      }
      
      Set<Int> points_on_facet(const Vector<Scalar>& hyp){
         Set<Int> result;
         for(Int i=0; i<points.rows(); i++){
            if(hyp*points.row(i) == 0) result += i;
         }
         return result;
      }

      Vector<Scalar> find_facet_hyperplane() {
         logger.log("find_facet_hyperplane", 2);
         Vector<Scalar> result(positive_hyperplane), next;
         Set<Int> facetPoints = points_on_facet(result);
         ListMatrix<Vector<Scalar>> frFacetPoints(0, points.cols());
         for(const auto& i : facetPoints){
            if(rank(frFacetPoints) < rank(frFacetPoints / points.row(i))){
               frFacetPoints /= points.row(i);
            }
         }
         ListMatrix<Vector<Scalar>> basis(0, points.cols());
         Matrix<Scalar> facet_affine_hull(orth_affine_hull);
         if(frFacetPoints.rows() > 0){
            facet_affine_hull /= null_space(frFacetPoints);
         }
         basis /= positive_hyperplane;
         bool has_next = true;
         while(has_next){
            // cout << "current result: " << result << endl;
            // cout << "points on facet: " << facetPoints << endl;
            next = find_next_basis_vector(result, facet_affine_hull);
            // cout << "next: " << next << endl;
            Scalar tilt_factor = find_tilting_factor(result, next);
            result += tilt_factor * next;
            basis /= result;
            // cout << "tilt factor: " << tilt_factor << endl;
            // cout << "new result: " << result << endl;
            Set<Int> comp(sequence(0,points.rows()) - facetPoints);
            for(const auto& i : comp){
               if(result * points.row(i) == 0){
                  facetPoints += i;
                  if(rank(frFacetPoints) < rank(frFacetPoints / points.row(i))){
                     frFacetPoints /= points.row(i);
                  }
               }
            }
            facet_affine_hull = orth_affine_hull / frFacetPoints;
            if(rank(frFacetPoints) == dim-1) has_next = false;
         }

         
         logger.log("find_facet_hyperplane done", 2);
         ls.lift_facet(result);
         orthogonalize(result, orth_affine_hull);
         normalize_facet_vector(result);
         logger.log(result, 2);
         return result;
      }

      void dualize_facet(const Vector<Scalar>& facet){
         logger.log("dualize_facet", 1);
         logger.log(facet, 1);
         // cout << orth_affine_hull << endl;
         Set<Int> facetPointsIndices(points_on_facet(facet));
         Matrix<Scalar> facetPoints(points.minor(facetPointsIndices, All));
         DCCH<Scalar> inner(facetPoints, positive_hyperplane, threshold, logger, ls, known + facets);
         Matrix<Scalar> facetFacets = inner.dualize();
         ls.push(points, orth_affine_hull, facet);
         for(const auto& ff : rows(facetFacets)){
            Vector<Scalar> tmp(ff);
            queue.push(ls.lift_facet(tmp));
         }
         ls.pop();
         logger.log("dualize_facet done", 1);
      }

      Vector<Scalar> lift_facet(const Vector<Scalar>& facetFacet, const Vector<Scalar>& facet){
         // logger.log("lift_facet");
         // logger.log(facetFacet);
         Vector<Scalar> result(facetFacet);
         // logger.log(result);
         Scalar factor(0);
         bool found = false;
         for(const auto& pt : rows(points)){
            Scalar fval = pt*facet;
            if(fval != 0){
               Scalar check = -(result*pt) / fval;
               if(!found){
                  factor = check;
                  found = true;
               }
               else if(check > factor){
                  factor = check;
               }
            }
         }
         // logger.log("Factor is: ");
         // cout << factor << endl;
         // logger.log("Facet is: ");
         // logger.log(facet);
         result += factor*facet;
         for(const auto& row : rows(orth_affine_hull)){
            // logger.log("Subtract row");
            // logger.log(Vector<Scalar>(row));
            Scalar rowval = row*row;
            if(rowval != 0){
               result -= ((result*row)/rowval) * row;
            }
            // logger.log(result);
         }
         normalize_facet_vector(result);
         // logger.log(result);
         // logger.log("lift_facet done");
         return result;
      }

      void dualize_recursion(){
         logger.log("dualize_recursion", 0);
         logger.descend();
         while(queue.size() != 0){
            logger.log("Queue size: " + std::to_string(queue.size()) + " (" + std::to_string(points.rows()) + ")", 0);
            Vector<Scalar> current = queue.front();
            queue.pop();
            if(!facets.contains(current) && !known.contains(current)){
               dualize_facet(current);
               facets += current;
               known += current;
            // } else {
            //    logger.log("Facet known");
            //    logger.log(current);
            }
            if(known.contains(current) && !facets.contains(current)){
               logger.log("This facet was known.", 3);
            }
         }
         logger.backtrack();
         logger.log("dualize_recursion done", 0);
      }

      void init() {
         logger.log("Building new dcch (" + std::to_string(points.rows()) + ")", 0);
         affine_hull = Set<Int>();
         std::string ahs("affine hull: ");
         for(Int i=0; i<points.rows(); i++){
            if(rank(points.minor(affine_hull, All)) < rank(points.minor(affine_hull + i, All))){
               affine_hull += i;
               ahs += std::to_string(i) + " ";
            }
         }
         orth_affine_hull = null_space(points.minor(affine_hull, All));
         logger.log(ahs, 1);
         for(Int i=0; i<orth_affine_hull.rows(); i++){
            const Vector<Scalar>& top(orth_affine_hull.row(i));
            Scalar toptop = top*top;
            for(Int j=i+1; j<orth_affine_hull.rows(); j++){
               Scalar bottomtop = orth_affine_hull.row(j)*top;
               orth_affine_hull.row(j) -= (bottomtop / toptop)*top;
            }
         }
         // cout << "orth_affine_hull:" << endl << orth_affine_hull << endl;
      }

   public:
      DCCH(const Matrix<Scalar>& pts, const Vector<Scalar>& pos_hyp, Int t, DCCH_Logger& l, LiftingStack<Scalar>& ls_in): 
         points(pts),
         positive_hyperplane(pos_hyp),
         threshold(t),
         logger(l),
         ls(ls_in)
         {
            init();
         }
      
      DCCH(const Matrix<Scalar>& pts, const Vector<Scalar>& pos_hyp, Int t, DCCH_Logger& l, LiftingStack<Scalar>& ls_in, const Set<Vector<Scalar>>& known_in): 
         points(pts),
         positive_hyperplane(pos_hyp),
         threshold(t),
         logger(l),
         ls(ls_in),
         known(known_in)
         {
            // logger.log("known facets: " + std::to_string(known.size()));
            init();
         }

      Matrix<Scalar> dualize(){
         logger.log("Dualizing", 0);
         Matrix<Scalar> result;
         if(points.rows() < threshold){
            logger.log("Below threshold " + std::to_string(threshold) + ", using enumerate_facets", 0);
            const auto facetData = enumerate_facets(points, zero_matrix<Scalar>(0, points.cols()), true);
            result = facetData.first;
         } else {
            logger.log("Above threshold " + std::to_string(threshold), 0);
            logger.descend();
            logger.log("Descended", 1);
            // cout << points << endl;
            // Maybe not compute this?
            dim = rank(points);
            logger.log("Dim is: " + std::to_string(dim), 0);
            Vector<Scalar> facet_vector;
            facet_vector = find_facet_hyperplane();
            queue.push(facet_vector);
            dualize_recursion();
            if(facets.size() > 0){
               result = Matrix<Scalar>(facets);
            } else {
               result = Matrix<Scalar>(0, points.cols());
            }
            logger.backtrack();
         }
         return result;
      }
};

} // namespace

template<typename Scalar>
Matrix<Scalar> dcch(const Matrix<Scalar>& P, const Vector<Scalar>& h, Int t, OptionSet options){
   const Int verbose = options["verbose"];
   // cout << "Log level is : " << verbose << endl;
   DCCH_Logger logger(verbose);
   LiftingStack<Scalar> ls;
   // cout << "Points:" << endl << P << endl;
   Matrix<Scalar> points(P);
   DCCH<Scalar> Dualizer(points, h, t, logger, ls);
   Matrix<Scalar> facets = Dualizer.dualize();
   // cout << "Facets are:" << endl << facets << endl;
   return facets;
}


UserFunctionTemplate4perl("# no doc yet",
                     "dcch<Scalar>(Matrix<type_upgrade<Scalar>>, Vector<type_upgrade<Scalar>>, Int; {verbose=>0})");

} // namespace polytope
} // namespace polymake
