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
Matrix<int> finish_set_monte_carlo(int dim, int tries, perl::OptionSet options);
Matrix<int> finish_set_n_moves(int dim, int tries, bool add, int n, int d);
Matrix<int> finish_given_set(Matrix<int> selected, int tries, perl::OptionSet options);
int compare_lex(Vector<int> a, Vector<int> b);


int compare_lex(Vector<int> a, Vector<int> b){
   Vector<int> diff = a - b;
   int i;
   for(i = 0; i<diff.dim(); i++){
      if(diff[i] < 0){
         return 1;
      }
      if(diff[i] > 0){
         return -1;
      }
   }
   return 0;
}

/////////////////////////////////////////////////
// Classes
//
class SetGame{
private:
   int dim;
   Map<Vector<int>, int> occupation;
   Matrix<int> possibles, selected, order;
   Vector<int> run_monte_carlo_try(Matrix<int> selected);

public:
   SetGame(int n, Matrix<int> o){
      dim = n;
      order = o;
      init();
   }

   SetGame(Matrix<int> preselected, Matrix<int> o){
      dim = preselected.cols();
      order = o;
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
      occupation = Map<Vector<int>, int>();
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

   int get_min_lines(){
      int result = 10000000;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            result = std::min(result, pt->second);
         }
      }
      return result;
   }

   int get_deviants(int minlines){
      int result = 0;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            if(pt->second != minlines){
               result++;
            }
         }
      }
      return result;
   }

   Vector<int> evaluate_board(){
      // return selected.rows();
      Vector<int> result(3);
      result[0] = selected.rows();
      result[1] = get_min_lines();
      result[2] = get_deviants(result[1]);
      return result;
   }

   int eval_n_random_moves(int n, int desired){
      int i = 0, presize = selected.rows(), result;
      while( (i < n) && (possibles.rows() > 0)){
         random_move();
         i++;
      }
      result = -get_deviants(desired);
      while(selected.rows() > presize){
         undo_last_no_update();
      }
      update_possibles();
      return result;
   }

   Vector<int> monte_carlo_try() {
      int presize = selected.rows();
      bool r = finish_monte_carlo();
      if(!r){ cout << "Something went wrong while playing!" << endl; }
      Vector<int> result = evaluate_board();
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
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second == 0){
            possibles /= pt->first;
         }
      }
   }

   void print(){
      cout << "Possibles: " << possibles.rows() << " Selected: " << selected.rows();
      if(possibles.rows() == 0){
         int minlines = get_min_lines();
         int deviants = get_deviants(minlines);
         cout << " minlines: " << minlines << " deviants: " << deviants;
      }
      cout << endl;
   }

   Vector<int> monte_carlo_evaluate(Vector<int> pt, int tries){
      move(pt);
      int i;
      Vector<int> result = monte_carlo_try();
      Vector<int> test;
      for(i=0; i<tries; i++){
         // cout << "Try: " << i << endl;
         // print();
         test = monte_carlo_try();
         if(compare_lex(order * result, order * test) == 1){
            result = test;
         }
         // result += monte_carlo_try();
         // print();
      }
      undo_last();
      return result;
   }

   Vector<int> find_best_move_monte_carlo(int tries){
      Matrix<int> pts(possibles);
      Matrix<int> currentBest(0,dim);
      Vector<int> currentBestVal = monte_carlo_try(), test;
      int i=0, total = pts.rows();
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(pts)); !pt.at_end(); pt++) {
         cout << "Point no. " << i << " of " << total;
         test = monte_carlo_evaluate(*pt, tries);
         cout << " gives " << test[0] << " - ";
         if(compare_lex(order*currentBestVal, order*test) == 1){
            cout << " WINNER!  ";
            currentBestVal = test;
            currentBest = Matrix<int>(0,dim);
            currentBest /= *pt;
         }
         if(compare_lex(order*currentBestVal, order*test) == 0){
            currentBest /= *pt;
         }
         cout << std::flush;
         i++;
      }
      cout << "Winning value is: " << currentBestVal << endl;
      i = rand() % currentBest.rows();
      cout << "Selecting index " << i << endl;
      return currentBest[i];
   }
   
   int n_moves_evaluate(Vector<int> pt, int tries, int n, int d){
      move(pt);
      int i;
      int result = INT_MIN;
      for(i=0; i<tries; i++){
         // cout << "Try: " << i << endl;
         // print();
         result = std::max(result, eval_n_random_moves(n, d));
         // result += monte_carlo_try();
         // print();
      }
      undo_last();
      return result;
   }
   
   Vector<int> find_best_move_n_moves(int tries, int n, int d){
      Matrix<int> pts(possibles);
      Matrix<int> currentBest(0,dim);
      int currentBestVal = INT_MIN, test, i=0, total = pts.rows();
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(pts)); !pt.at_end(); pt++) {
         cout << "Point no. " << i << " of " << total;
         test = n_moves_evaluate(*pt, tries, n, d);
         cout << " gives " << test << " - ";
         if(test > currentBestVal){
            cout << " WINNER!  ";
            currentBestVal = test;
            currentBest = Matrix<int>(0,dim);
            currentBest /= *pt;
         }
         if(test == currentBestVal){
            currentBest /= *pt;
         }
         cout << std::flush;
         i++;
      }
      cout << "Winning value is: " << currentBestVal << endl;
      i = rand() % currentBest.rows();
      cout << "Selecting index " << i << endl;
      return currentBest[i];
   }

   void finish_game_monte_carlo(int tries, bool add){
      Vector<int> next;
      int additional = 0;
      while(possibles.rows() > 0){
         print();
         next = find_best_move_monte_carlo(tries + additional);
         // print();
         additional += add ? 1 : 0;
         move(next);
      }
      print();
   }
   
   void finish_game_n_moves(int tries, bool add, int n, int d){
      Vector<int> next;
      int additional = 0;
      while(possibles.rows() > 0){
         print();
         next = find_best_move_n_moves(tries + additional, n, d);
         // print();
         additional += add ? 1 : 0;
         move(next);
      }
      print();
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
   
Vector<int> run_monte_carlo_try(Matrix<int> selected){
   SetGame set = SetGame(selected, unit_matrix<int>(3));
   set.print();
   Vector<int> result = set.monte_carlo_try();
   set.print();
   return result;
}

Matrix<int> finish_given_set(Matrix<int> selected, int tries, perl::OptionSet options){
   bool add = options["add"];
   Matrix<int> order = options["order"];
   SetGame set = SetGame(selected, order);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<int> finish_set_monte_carlo(int dim, int tries, perl::OptionSet options){
   int i;
   bool add = options["add"];
   Matrix<int> order = options["order"];
   SetGame set = SetGame(dim, order);
   Vector<int> v = zero_vector<int>(dim);
   set.move_no_update(v);
   for(i=0; i<dim; i++){
      v = unit_vector<int>(dim, i);
      set.move_no_update(v);
   }
   set.update_possibles();
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<int> finish_set_n_moves(int dim, int tries, bool add, int n, int d){
   SetGame set = SetGame(dim, unit_matrix<int>(3));
   int i;
   Vector<int> v = zero_vector<int>(dim);
   set.move_no_update(v);
   for(i=0; i<dim; i++){
      v = unit_vector<int>(dim, i);
      set.move_no_update(v);
   }
   set.update_possibles();
   set.finish_game_n_moves(tries, add, n, d);
   return set.get_selected();
}


} // namespace

Function4perl(&get_init_matrix, "get_init_matrix");

Function4perl(&get_third_on_line, "get_third_on_line");

Function4perl(&run_monte_carlo_try, "run_monte_carlo_try");

Function4perl(&finish_set_monte_carlo, "finish_set_monte_carlo( $ , $ , {options=>1, order=>unit_matrix<Int>(3)})");

Function4perl(&finish_set_n_moves, "finish_set_n_moves");

Function4perl(&finish_given_set, "finish_given_set( Matrix<Int> , $ , {options=>1, order=>unit_matrix<Int>(3)})");

} // namespace polymake
} // namespace polytope
