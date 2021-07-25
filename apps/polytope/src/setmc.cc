/* Copyright (c) 1997-2017
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
#include "map"
#include "polymake/linalg.h"
#include "polymake/vector"
#include "polymake/PowerSet.h"
#include "polymake/common/lattice_tools.h"


namespace polymake {
namespace polytope {
namespace {


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunused-private-field"
#pragma clang diagnostic ignored "-Wreorder"
#endif



Int EC = 0;

Matrix<Int> get_init_matrix(Int dim);
Vector<Int> get_third_on_line(Vector<Int> a, Vector<Int> b);
Matrix<Int> run_monte_carlo_try(Matrix<Int> selected);
Matrix<Int> run_monte_carlo_try_symmetric(Matrix<Int> selected, Matrix<Int> sym, Vector<Int> shift);
Matrix<Int> finish_set(Int dim, Int tries, OptionSet options);
Matrix<Int> finish_inverse_set(Int dim);
Matrix<Int> finish_given_set(Matrix<Int> selected, Int tries, OptionSet options);
Matrix<Int> finish_given_set_symmetric(Matrix<Int> selected, Int tries, Matrix<Int> A, Vector<Int> b, OptionSet options);
Int compare_lex(Vector<Int> a, Vector<Int> b);
Int scalp(Vector<Int> a, Vector<Int> b);
Vector<Int> tof3(Vector<Int> a);
Matrix<Int> points_on_hyperplane(Matrix<Int> given, Vector<Int> normal, Int val);
ListReturn find_worst_hyperplane(Matrix<Int> given);

void vectorMod3(Vector<Int>& b){
   for(Int i=0; i<b.dim(); i++){
      b[i] %= 3;
      b[i] += 3;
      b[i] %= 3;
   }
}

Vector<Int> get_third_on_line(Vector<Int> a, Vector<Int> b){
   Vector<Int> result = -a-b;
   vectorMod3(result);
   return result;
}


Int compare_lex(Vector<Int> a, Vector<Int> b){
   Vector<Int> diff = a - b;
   Int i;
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
// Inverse Set
//
class InverseSetGame{
private:
   Int dim;
   Map<Vector<Int>, Int> linesThroughPoint;
   Matrix<Int> selected, good_ones;

public:
   InverseSetGame(Int n){
      dim = n;
      init();
   }

   void init() {
      Int lines;
      good_ones = get_init_matrix(dim);
      lines = (good_ones.rows() - 1)/2;
      selected = Matrix<Int>(0, dim);
      linesThroughPoint = Map<Vector<Int>, Int>();
      for(const auto& g : rows(good_ones)) {
         linesThroughPoint[g] = lines;
      }
   }

   bool move(Vector<Int> v){
      bool result = move_no_update(v);
      if(!result){ return result;}
      update_good_ones();
      return result;
   }

   bool move_no_update(Vector<Int> v){
      Vector<Int> third;
      if(linesThroughPoint[v] == -2){
         return false;
      }
      for(const auto& pt : linesThroughPoint) {
         if(pt.second >= 0){
            linesThroughPoint[pt.first]--;
         }
      }
      for(const auto& g : rows(selected)) {
         third = get_third_on_line(g, v);
         if(linesThroughPoint[third] != -2){
            linesThroughPoint[third]++;
         }
      }
      linesThroughPoint[v] = -2;
      selected /= v;
      return true;
   }

   Int find_max_val(){
      Int result = 0;
      for(const auto& pt : linesThroughPoint) {
         if(pt.second > result){
            result = pt.second;
         }
      }
      return result;
   }

   Matrix<Int> find_best_moves(){
      Int currentBestVal = find_max_val();
      Matrix<Int> result(0,dim);
      for(const auto& pt : linesThroughPoint) {
         if(pt.second == currentBestVal){
            result /= pt.first;
         }
      }
      return result;
   }

   bool random_best_move(){
      Matrix<Int> possibles = find_best_moves();
      Int index = rand() % possibles.rows();
      return move(possibles[index]);
   }
   
   bool random_move(){
      Int index = rand() % good_ones.rows();
      return move(good_ones[index]);
   }

   Matrix<Int> finish_randomly(){
      while(find_max_val() > 0){
         random_move();
      }
      return selected;
   }

   void update_good_ones(){
      good_ones = Matrix<Int>(0, dim);
      for(const auto& pt : linesThroughPoint) {
         if(pt.second > 0){
            good_ones /= pt.first;
         }
      }  
   }
};

/////////////////////////////////////////////////
// Set
//
class SetGame{
private:
   Int dim, lookahead, targetdev, symmc;
   Map<Vector<Int>, Int> occupation;
   Matrix<Int> possibles, selected, order;
   bool symmetric;
   Matrix<Int> symmMat;
   Vector<Int> symmShift;

public:
   SetGame(Int n, Matrix<Int> o){
      dim = n;
      init(o);
   }

   SetGame(Matrix<Int> preselected, Matrix<Int> o){
      dim = preselected.cols();
      init(o);
      for(const auto& s : rows(preselected)) {
         move_no_update(s);
      }
      update_possibles();
   }

   void set_symmetry(Matrix<Int> sM, Vector<Int> b){
      symmetric = true;
      symmMat = sM;
      symmShift = b;
   }

   void set_target_deviation(Int t){
      targetdev = t;
   }

   Matrix<Int> get_selected(){
      return Matrix<Int>(selected);
   }
   
   void init(Matrix<Int> o) {
      symmetric = false;
      symmc = 0;
      lookahead = INT_MAX;
      targetdev = -1;
      order = o;
      possibles = get_init_matrix(dim);
      selected = Matrix<Int>(0, dim);
      occupation = Map<Vector<Int>, Int>();
      for(const auto& g : rows(possibles)) {
         occupation[g] = 0;
      }
   }

   bool move_no_update(Vector<Int> v){
      Vector<Int> third, partner;
      if(symmetric){
         partner = tof3(symmMat*v + symmShift);
      }
      if(occupation[v] != 0){
         // cout << "Returning false." << endl;
         return false;
      }
      occupation[v] = -1;
      for(const auto& s : rows(selected)) {
         third = get_third_on_line(s, v);
         // if(occupation[third] == -1){
         //    cout << endl << "Something went wrong! " << third << " - " << *s << " - " << v << endl;
         // }
         occupation[third]++;
      }
      selected /= v;
      if(symmetric && (symmc < 1000) && (occupation[partner] != -1)){
         //cout << v << " - " << partner << endl;
         if(partner != v){
            bool s = move_no_update(partner);
            if(!s){
               // cout << "Symmetric move was not possible." << endl;
               symmc++;
               undo_last_no_update();
               return s;
            }
         }
      } else if (symmetric){
         // cout << "Resetting symmc." << endl;
         symmc = 0;
      }
      return true;
   }

   Int get_min_lines(){
      Int result = INT_MAX;
      for(const auto& pt : occupation) {
         if(pt.second != -1){
            result = std::min(result, pt.second);
         }
      }
      return result;
   }
   
   Int get_max_lines(){
      Int result = INT_MIN;
      for(const auto& pt : occupation) {
         if(pt.second != -1){
            result = std::max(result, pt.second);
         }
      }
      return result;
   }
   
   Int get_total_lines(){
      Int result = 0;
      for(const auto& pt : occupation) {
         if(pt.second != -1){
            result += pt.second;
         }
      }
      return result;
   }
   

   Int get_deviants(Int minlines){
      Int result = 0;
      for(const auto& pt : occupation) {
         if(pt.second != -1){
            if(pt.second != minlines){
               result += pt.second == 0 ? 0 : 1;
            }
         }
      }
      return result;
   }

   Vector<Int> evaluate_board(){
      // return selected.rows();
      Vector<Int> result(6);
      result[0] = selected.rows();
      result[1] = get_min_lines();
      if(targetdev == -1){
         result[2] = get_deviants(result[1]);
      } else {
         result[2] = get_deviants(targetdev);
      }
      result[3] = get_max_lines();
      result[4] = get_total_lines();
      result[5] = possibles.rows();
      return result;
   }
   
   void print(){
      Vector<Int> eval = evaluate_board();
      cout << "Current evaluation:" << endl;
      cout << "Selected(0): " << eval[0] << " ";
      cout << "Possibles(5): " << eval[5] << " ";
      cout << "Deviants(2): " << eval[2] << " ";
      cout << "Minlines(1): " << eval[1] << " ";
      cout << "Maxlines(3): " << eval[3] << " ";
      cout << "TotalLines(4): " << eval[4] << " ";
      cout << endl;
   }


   Vector<Int> monte_carlo_try_no_undo() {
      bool r = finish_monte_carlo();
      if(!r){ cout << "Something went wrong while playing!" << endl; }
      return evaluate_board();
   }

   Vector<Int> monte_carlo_try(){
      Int presize = selected.rows();
      Vector<Int> result = monte_carlo_try_no_undo();
      while(selected.rows() > presize){
         undo_last_no_update();
      }
      update_possibles();
      return result;
   }

   void set_lookahead(Int n){
      lookahead = n;
   }

   bool finish_monte_carlo(){
      Int i = 0;
      while((possibles.rows() > 0) && (i<lookahead)){
         bool r = random_move();
         if(!r && (EC > 10000)){ 
            EC = 0;
            return r; 
         } else if(!r) {
            EC++;
         }
         i++;
      }
      return true;
   }

   bool random_move(){
      Int i = rand() % possibles.rows();
      return move(possibles[i]);
   }

   bool move(Vector<Int> v){
      bool result = move_no_update(v);
      update_possibles();
      return result;
   }

   bool undo_last_no_update(){
      Int index = selected.rows() - 1;
      if(index == -1){
         return false;
      }
      Vector<Int> last = selected[index], third;
      selected = selected.minor(~scalar2set(index), All);
      occupation[last] = 0;
      for(const auto& s : rows(selected)) {
         third = get_third_on_line(s, last);
         // if(occupation[third] <= 0){
         //    cout << endl << "Something went wrong undo!" << endl;
         // }
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
      possibles = Matrix<Int>(0, dim);
      for(const auto& pt : occupation) {
         if(pt.second == 0){
            possibles /= pt.first;
         }
      }
   }


   Vector<Int> monte_carlo_evaluate(Vector<Int> pt, Int tries){
      Int i, presize=selected.rows();
      move(pt);
      Vector<Int> result = monte_carlo_try();
      Vector<Int> test;
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
      while(selected.rows() > presize){
         undo_last();
      }
      return result;
   }

   Vector<Int> find_best_move_monte_carlo(Int tries){
      Matrix<Int> pts(possibles);
      // if(pts.rows() == 0){
      //    cout << endl << "I went here, even though there are no possibles." << endl;
      // }
      Matrix<Int> currentBest(0,dim);
      currentBest /= possibles[0];
      Vector<Int> currentBestVal = monte_carlo_evaluate(possibles[0], tries), test;
      Int i=0, total = pts.rows();
      for(const auto& pt : rows(pts)) {
         cout << "Point no. " << i << " of " << total;
         test = monte_carlo_evaluate(pt, tries);
         cout << " gives " << test[0] << " - ";
         // if(occupation[*pt] != 0){
         //    cout << endl << "Something went wrong! " << occupation[*pt] << endl;
         // }
         if(compare_lex(order*currentBestVal, order*test) == 1){
            cout << " WINNER!  ";
            currentBestVal = test;
            currentBest = Matrix<Int>(0,dim);
            currentBest /= pt;
         }
         if(compare_lex(order*currentBestVal, order*test) == 0){
            currentBest /= pt;
         }
         cout << std::flush;
         i++;
      }
      cout << "Winning value is: " << currentBestVal << " : " << currentBest.rows() << endl;
      i = currentBest.rows() == 1 ? 0 : rand() % currentBest.rows();
      cout << "Selecting index " << i << " : " << currentBest[i] << endl;
      return currentBest[i];
   }
   

   void finish_game_monte_carlo(Int tries, bool add){
      Vector<Int> next;
      Int additional = 0;
      while(possibles.rows() > 0){
         print();
         next = find_best_move_monte_carlo(tries + additional);
         // print();
         additional += add ? 1 : 0;
         cout << "Doing next move." << endl;
         move(next);
      }
      print();
   }
   
};

/////////////////////////////////////////////////
// ReverseSet
//
class ReverseSet{
private:
   Int dim;
   Map<Vector<Int>, Int> pointIndices;
   Matrix<Int> points, lines;

public:
   ReverseSet(Int d){
      dim = d;
      points = get_init_matrix(d);
      // cout << "Points ok." << endl;
      lines(0,3);
      Vector<Int> third(3), line(3);
      Int i = 0, j, thirdIndex;
      Int n = points.rows();
      for(const auto& pt : rows(points)) {
         pointIndices[pt] = i;
         i++;
      }
      // cout << "PointIndices ok." << endl;
      for(i=0; i<n; i++){
         for(j=i+1; j<n; j++){
            // cout << "Computing third. " << points[i] << "," << points[j] << endl;
            third = tof3(get_third_on_line(points[i], points[j]));
            // cout << "Third is: " << third << " - " << points[i] << "," << points[j] << endl;
            thirdIndex = pointIndices[third];
            if(j<thirdIndex){
               line[0] = i;
               line[1] = j;
               line[2] = thirdIndex;
               // cout << "Adding line. " << line << endl;
               lines = lines/Vector<Int>(line);
            }
         }
      }
   }

   ReverseSet(Matrix<Int> pts, Map<Vector<Int>, Int> pI, Matrix<Int> l){
      points = pts;
      pointIndices = pI;
      lines = l;
   }

   Matrix<Int> get_points(){
      return points;
   }

   Matrix<Int> get_lines(){
      return lines;
   }

   Map<Vector<Int>, Int> get_pt_index_map(){
      return pointIndices;
   }

   void remove_point(const Vector<Int> pt){
      Int pointIndex = pointIndices[pt];
      points = points.minor(~scalar2set(find_index_of_pt(pt)), All);
      Matrix<Int> keepLines(0,3);
      for(const auto& line : rows(lines)) {
         if(!line_contains(line, pointIndex)){
            keepLines /= line;
         }
      }
      lines = keepLines;
   }

   bool line_contains(const Vector<Int> line, Int index){
      if(line[0] == index){
         return true;
      }
      if(line[1] == index){
         return true;
      }
      if(line[2] == index){
         return true;
      }
      return false;
   }

   Int find_index_of_pt(const Vector<Int> p){
      Int i = 0;
      for(const auto& pt : rows(points)) {
         if(p == pt){
            return i;
         }
         i++;
      }
      return -1;
   }
};

//////////////////////////////////////////////////
// Methods
//


Int scalp(Vector<Int> a, Vector<Int> b){
   Int result = 0, i;
   for(i=0; i<a.dim(); i++){
      result += a[i] * b[i];
   }
   return result % 3;
}

Matrix<Int> points_on_hyperplane(Matrix<Int> given, Vector<Int> normal, Int val){
   Int dim = given.cols();
   Matrix<Int> result(0,dim);
   for(const auto& pt : rows(given)) {
      Int test = scalp(pt, normal);
      if(test == val){
         result /= pt;
      }
   }
   return result;
}


ListReturn get_init_reverse_set(Int d){
   ReverseSet rs(d);
   ListReturn result;
   result << rs.get_points();
   result << rs.get_pt_index_map();
   result << rs.get_lines();
   return result;
}

ListReturn remove_pt_from_reverse_set(Matrix<Int> pts, Map<Vector<Int>, Int> pI, Matrix<Int> lines, Vector<Int> pt){
   ReverseSet rs(pts, pI, lines);
   rs.remove_point(pt);
   ListReturn result;
   result << rs.get_points();
   result << rs.get_pt_index_map();
   result << rs.get_lines();
   return result;
}

ListReturn find_worst_hyperplane(Matrix<Int> given){
   Int dim = given.cols(), i, min = given.rows() + 1, val = 0;
   Matrix<Int> hyperplanes = get_init_matrix(dim), A;
   hyperplanes = hyperplanes.minor(~scalar2set(0),All);
   Vector<Int> hp;
   ListReturn result;
   for(const auto& n : rows(hyperplanes)) {
      for(i=0; i<3; i++){
         A = points_on_hyperplane(given, n, i);
         cout << "Vector: " << n << " val: " << i << " gives " << A.rows() << endl;
         if(A.rows() < min){
            min = A.rows();
            hp = n;
            val = i;
         }
      }
   }
   result << hp << val << points_on_hyperplane(given, hp, val);
   return result;
}

Matrix<Int> get_init_matrix(Int dim){
   if(dim == 1){
      Matrix<Int> result(3,1);
      result(0,0) = 0;
      result(1,0) = 1;
      result(2,0) = 2;
      return result;
   }
   Matrix<Int> preresult = get_init_matrix(dim-1);
   Matrix<Int> result = preresult | zero_vector<Int>(preresult.rows());
   result /= preresult | ones_vector<Int>(preresult.rows());
   result /= preresult | (2 * ones_vector<Int>(preresult.rows()));
   return result;
}

   
Matrix<Int> run_monte_carlo_try(Matrix<Int> selected){
   SetGame set = SetGame(selected, unit_matrix<Int>(6));
   // set.print();
   Vector<Int> result = set.monte_carlo_try_no_undo();
   if(result[0] == 15){
      cout << "Blip!" << endl;
   }
   // set.print();
   return set.get_selected();
}

Matrix<Int> run_monte_carlo_try_symmetric(Matrix<Int> selected, Matrix<Int> sym, Vector<Int> shift){
   SetGame set = SetGame(selected, unit_matrix<Int>(6));
   set.set_symmetry(sym, shift);
   // set.print();
   Vector<Int> result = set.monte_carlo_try_no_undo();
   if(result[0] == 15){
      cout << "Blip!" << endl;
   }
   // set.print();
   return set.get_selected();
}

Matrix<Int> finish_given_set(Matrix<Int> selected, Int tries, OptionSet options){
   Int la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<Int> order = options["order"];
   SetGame set = SetGame(selected, order);
   set.set_lookahead(la);
   set.set_target_deviation(t);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<Int> finish_given_set_symmetric(Matrix<Int> selected, Int tries, Matrix<Int> A, Vector<Int> b, OptionSet options){
   Int la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<Int> order = options["order"];
   SetGame set = SetGame(selected, order);
   set.set_symmetry(A,b);
   set.set_lookahead(la);
   set.set_target_deviation(t);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<Int> finish_set(Int dim, Int tries, OptionSet options){
   Int i, la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<Int> order = options["order"];
   SetGame set = SetGame(dim, order);
   set.set_lookahead(la);
   set.set_target_deviation(t);
   Vector<Int> v = zero_vector<Int>(dim);
   set.move_no_update(v);
   for(i=0; i<dim; i++){
      v = unit_vector<Int>(dim, i);
      set.move_no_update(v);
   }
   set.update_possibles();
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Vector<Int> tof3(Vector<Int> a){
   Int i;
   Vector<Int> result(a);
   for(i=0; i<result.dim(); i++){
      result[i] = (result[i]+300) % 3;
   }
   return result;
}

Matrix<Int> finish_inverse_set(Int dim){
   InverseSetGame IS(dim);
   return IS.finish_randomly();
}


#if defined(__clang__)
#pragma clang diagnostic pop
#endif
} // namespace

Function4perl(&get_init_reverse_set, "get_init_reverse_set");

Function4perl(&remove_pt_from_reverse_set, "remove_pt_from_reverse_set");

Function4perl(&get_init_matrix, "get_init_matrix");

Function4perl(&get_third_on_line, "get_third_on_line");

Function4perl(&run_monte_carlo_try, "run_monte_carlo_try");

Function4perl(&run_monte_carlo_try_symmetric, "run_monte_carlo_try_symmetric");

Function4perl(&find_worst_hyperplane, "find_worst_hyperplane");

Function4perl(&finish_inverse_set, "finish_inverse_set");

Function4perl(&finish_set, "finish_set( $ , $ , {add=>1, order=>(unit_matrix<Int>(6)), targetdev=>-1, lookahead=>10000})");

Function4perl(&finish_given_set, "finish_given_set( Matrix<Int> , $ , {add=>1, order=>(unit_matrix<Int>(6)), targetdev=>-1, lookahead=>10000})");

Function4perl(&finish_given_set_symmetric, "finish_given_set_symmetric( Matrix<Int> , $ , Matrix<Int>, Vector<Int>, {add=>1, order=>(unit_matrix<Int>(6)), targetdev=>-1, lookahead=>10000})");

} // namespace polymake
} // namespace polytope
