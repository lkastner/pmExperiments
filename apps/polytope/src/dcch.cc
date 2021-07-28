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
#include <random>



namespace polymake {
namespace polytope {
namespace {

class DCCH_Logger {
   private:
      Int depth;

   public:
      DCCH_Logger() : depth(0) {}

      void descend() {
         depth++;
      }

      void backtrack() {
         depth--;
      }

      void log(std::string text) {
         for(int i=0; i<depth; i++){
            cout << "| ";
         }
         cout << text << endl;
      }
};

template<typename Scalar>
class DCCH {
   private:
      const Matrix<Scalar>& points;
      const Vector<Scalar>& positive_hyperplane;
      const Int threshold;
      DCCH_Logger& logger;
      Int dim;

      bool is_separating(const Vector<Scalar>& h){
         bool negative=false, positive=false;
         for(const auto& p : points){
            if(p*h < 0){
               negative = true;
            }
            if(p*h > 0){
               positive = true;
            }
            if(positive && negative) return true;
         }
         return false;
      }

      Vector<Scalar> find_facet_hyperplane() {
         logger.log("find_facet_hyperplane");
         std::random_device rd; // obtain a random number from hardware
         std::mt19937 gen(rd()); // seed the generator
         std::uniform_int_distribution<> distr(0, points.rows()); // define the range

         Vector<Scalar> start(points[distr(gen)] - points[distr(gen)]);
         if(is_zero(start)){
            throw std::runtime_error("Identical pts");
         } else {
            
         }
         if(is_separating(start)){
            logger.log("Was separating");
            Matrix<Scalar> basis(complete_to_orth_basis(start));
         }
         return start;
      }

   public:
      DCCH(const Matrix<Scalar>& pts, const Vector<Scalar>& pos_hyp, Int t, DCCH_Logger& l): 
         points(pts),
         positive_hyperplane(pos_hyp),
         threshold(t),
         logger(l)
         {
            logger.log("Building new dcch");
         }

      void dualize(){
         logger.log("Dualizing");
         logger.descend();
         // Maybe not compute this?
         dim = rank(points);
         logger.log("Dim is: " + std::to_string(dim));
         Vector<Scalar> facet_vector = find_facet_hyperplane();

         logger.backtrack();
      }
};

} // namespace

template<typename Scalar>
void dcch(const Matrix<Scalar>& P, const Vector<Scalar>& h){
   DCCH_Logger logger;
   DCCH<Scalar> Dualizer(P, h, 100, logger);
   Dualizer.dualize();
}


UserFunctionTemplate4perl("# no doc yet",
                     "dcch<Scalar>(Matrix<type_upgrade<Scalar>>, Vector<type_upgrade<Scalar>>)");

} // namespace polymake
} // namespace polytope
