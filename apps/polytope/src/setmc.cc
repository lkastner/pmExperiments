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
#include "polymake/Map.h"
#include "polymake/vector"
#include "polymake/PowerSet.h"
#include "polymake/common/lattice_tools.h"


namespace polymake {
namespace polytope {
namespace {

Matrix<int> get_init_matrix(int dim);
Vector<int> get_third_on_line(Vector<int> a, Vector<int> b);
Matrix<int> finish_set(Matrix<int> selected, int tries, bool add);

/////////////////////////////////////////////////
// Classes
//
class SetGame{
private:
   int dim;
   Map<Vector<int>, Integer> occupation;
   Matrix<int> possibles, selected;
   int run_monte_carlo_try(Matrix<int> selected);

public:
   SetGame(int n){
      dim = n;
      init();
   }

   SetGame(Matrix<int> preselected){
      dim = preselected.cols();
      init();
      for(Entire<Rows<Matrix<int> > >::const_iterator s = entire(rows(preselected)); !s.at_end(); s++) {
         move_no_update(*s);
      }
      update_possibles();
   }

   Matrix<int> get_selected(){
      return Matrix<int>(selected);
   }
   
   void init() {
      possibles = get_init_matrix(dim);
      selected = Matrix<int>(0, dim);
      occupation = Map<Vector<int>, Integer>();
      for(Entire<Rows<Matrix<int> > >::const_iterator g = entire(rows(possibles)); !g.at_end(); g++) {
         occupation[*g] = 0;
      }
   }

   bool move_no_update(Vector<int> v){
      Vector<int> third;
      if(occupation[v] != 0){
         return false;
      }
      occupation[v] = -1;
      for(Entire<Rows<Matrix<int> > >::const_iterator s = entire(rows(selected)); !s.at_end(); s++) {
         third = get_third_on_line(*s, v);
         occupation[third]++;
      }
      selected /= v;
      return true;
   }

   int monte_carlo_try() {
      int presize = selected.rows();
      bool r = finish_monte_carlo();
      if(!r){ return r; }
      int result = selected.rows();
      while(selected.rows() > presize){
         undo_last_no_update();
      }
      update_possibles();
      return result;
   }

   bool finish_monte_carlo(){
      while(possibles.rows() > 0){
         bool r = random_move();
         if(!r){ return r; }
      }
      return true;
   }

   bool random_move(){
      int i = rand() % possibles.rows();
      return move(possibles[i]);
   }

   bool move(Vector<int> v){
      bool result = move_no_update(v);
      update_possibles();
      return result;
   }

   bool undo_last_no_update(){
      int index = selected.rows() - 1;
      if(index == -1){
         return false;
      }
      Vector<int> last = selected[index], third;
      selected = selected.minor(~scalar2set(index), All);
      occupation[last] = 0;
      for(Entire<Rows<Matrix<int> > >::const_iterator s = entire(rows(selected)); !s.at_end(); s++) {
         third = get_third_on_line(*s, last);
         occupation[third]--;
      }
      return true;
   }

   bool undo_last(){
      bool result = undo_last_no_update();
      update_possibles();
      return result;
   }

   void update_possibles(){
      possibles = Matrix<int>(0, dim);
      for(Entire<Map<Vector<int>, Integer> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second == 0){
            possibles /= pt->first;
         }
      }
   }

   void print(){
      cout << "Possibles: " << possibles.rows() << " Selected: " << selected.rows() << endl;
   }

   int monte_carlo_evaluate(Vector<int> pt, int tries){
      move(pt);
      int i;
      int result = 0;
      for(i=0; i<tries; i++){
         // cout << "Try: " << i << endl;
         // print();
         result = std::max(result, monte_carlo_try());
         // print();
      }
      undo_last();
      return result;
   }

   Vector<int> find_best_move(int tries){
      Matrix<int> pts(possibles);
      Vector<int> currentBest;
      int currentBestVal = 0, test, i=0, total = pts.rows();
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(pts)); !pt.at_end(); pt++) {
         cout << "Point no. " << i << " of " << total;
         test = monte_carlo_evaluate(*pt, tries);
         cout << " gives " << test << endl;
         if(test > currentBestVal){
            currentBest = *pt;
         }
         i++;
      }
      return currentBest;
   }

   void finish_game_monte_carlo(int tries, bool add){
      Vector<int> next;
      int additional = 0;
      while(possibles.rows() > 0){
         print();
         next = find_best_move(tries + additional);
         print();
         additional += add ? 1 : 0;
         move(next);
      }
   }
};

Matrix<int> get_init_matrix(int dim){
   if(dim == 1){
      Matrix<int> result(3,1);
      result(0,0) = 0;
      result(1,0) = 1;
      result(2,0) = 2;
      return result;
   }
   Matrix<int> preresult = get_init_matrix(dim-1);
   Matrix<int> result = preresult | zero_vector<int>(preresult.rows());
   result /= preresult | ones_vector<int>(preresult.rows());
   result /= preresult | (2 * ones_vector<int>(preresult.rows()));
   return result;
}

Vector<int> get_third_on_line(Vector<int> a, Vector<int> b){
   Vector<int> result = -a-b;
   int i;
   for(i=0; i<result.dim(); i++){
      result[i] += 6;
      result[i] %= 3;
   }
   return result;
}
   
int run_monte_carlo_try(Matrix<int> selected){
   SetGame set = SetGame(selected);
   set.print();
   int result = set.monte_carlo_try();
   set.print();
   return result;
}

Matrix<int> finish_set(Matrix<int> selected, int tries, bool add){
   SetGame set = SetGame(selected);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}


} // namespace

Function4perl(&get_init_matrix, "get_init_matrix");

Function4perl(&get_third_on_line, "get_third_on_line");

Function4perl(&run_monte_carlo_try, "run_monte_carlo_try");

Function4perl(&finish_set, "finish_set");

} // namespace polymake
} // namespace polytope
