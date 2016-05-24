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

int EC = 0;

Matrix<int> get_init_matrix(int dim);
Vector<int> get_third_on_line(Vector<int> a, Vector<int> b);
Matrix<int> run_monte_carlo_try(Matrix<int> selected);
Matrix<int> run_monte_carlo_try_symmetric(Matrix<int> selected, Matrix<int> sym, Vector<int> shift);
Matrix<int> finish_set(int dim, int tries, perl::OptionSet options);
Matrix<int> finish_inverse_set(int dim);
Matrix<int> finish_given_set(Matrix<int> selected, int tries, perl::OptionSet options);
Matrix<int> finish_given_set_symmetric(Matrix<int> selected, int tries, Matrix<int> A, Vector<int> b, perl::OptionSet options);
int compare_lex(Vector<int> a, Vector<int> b);
int scalp(Vector<int> a, Vector<int> b);
Vector<int> tof3(Vector<int> a);
Matrix<int> points_on_hyperplane(Matrix<int> given, Vector<int> normal, int val);
perl::ListReturn find_worst_hyperplane(Matrix<int> given);


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
// Inverse Set
//
class InverseSetGame{
private:
   int dim;
   Map<Vector<int>, int> linesThroughPoint;
   Matrix<int> selected, good_ones;

public:
   InverseSetGame(int n){
      dim = n;
      init();
   }

   void init() {
      int lines;
      good_ones = get_init_matrix(dim);
      lines = (good_ones.rows() - 1)/2;
      selected = Matrix<int>(0, dim);
      linesThroughPoint = Map<Vector<int>, int>();
      for(Entire<Rows<Matrix<int> > >::const_iterator g = entire(rows(good_ones)); !g.at_end(); g++) {
         linesThroughPoint[*g] = lines;
      }
   }

   bool move(Vector<int> v){
      bool result = move_no_update(v);
      if(!result){ return result;}
      update_good_ones();
      return result;
   }

   bool move_no_update(Vector<int> v){
      Vector<int> third;
      if(linesThroughPoint[v] == -2){
         return false;
      }
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(linesThroughPoint); !pt.at_end(); ++pt){
         if(pt->second >= 0){
            linesThroughPoint[pt->first]--;
         }
      }
      for(Entire<Rows<Matrix<int> > >::const_iterator g = entire(rows(selected)); !g.at_end(); g++) {
         third = get_third_on_line(*g, v);
         if(linesThroughPoint[third] != -2){
            linesThroughPoint[third]++;
         }
      }
      linesThroughPoint[v] = -2;
      selected /= v;
      return true;
   }

   int find_max_val(){
      int result = 0;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(linesThroughPoint); !pt.at_end(); ++pt){
         if(pt->second > result){
            result = pt->second;
         }
      }
      return result;
   }

   Matrix<int> find_best_moves(){
      int currentBestVal = find_max_val();
      Matrix<int> result(0,dim);
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(linesThroughPoint); !pt.at_end(); ++pt){
         if(pt->second == currentBestVal){
            result /= pt->first;
         }
      }
      return result;
   }

   bool random_best_move(){
      Matrix<int> possibles = find_best_moves();
      int index = rand() % possibles.rows();
      return move(possibles[index]);
   }
   
   bool random_move(){
      int index = rand() % good_ones.rows();
      return move(good_ones[index]);
   }

   Matrix<int> finish_randomly(){
      while(find_max_val() > 0){
         random_move();
      }
      return selected;
   }

   void update_good_ones(){
      good_ones = Matrix<int>(0, dim);
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(linesThroughPoint); !pt.at_end(); ++pt){
         if(pt->second > 0){
            good_ones /= pt->first;
         }
      }  
   }
};

/////////////////////////////////////////////////
// Set
//
class SetGame{
private:
   int dim, lookahead, targetdev, symmc;
   Map<Vector<int>, int> occupation;
   Matrix<int> possibles, selected, order;
   bool symmetric;
   Matrix<int> symmMat;
   Vector<int> symmShift;

public:
   SetGame(int n, Matrix<int> o){
      dim = n;
      init(o);
   }

   SetGame(Matrix<int> preselected, Matrix<int> o){
      dim = preselected.cols();
      init(o);
      for(Entire<Rows<Matrix<int> > >::const_iterator s = entire(rows(preselected)); !s.at_end(); s++) {
         move_no_update(*s);
      }
      update_possibles();
   }

   void set_symmetry(Matrix<int> sM, Vector<int> b){
      symmetric = true;
      symmMat = sM;
      symmShift = b;
   }

   void set_target_deviation(int t){
      targetdev = t;
   }

   Matrix<int> get_selected(){
      return Matrix<int>(selected);
   }
   
   void init(Matrix<int> o) {
      symmetric = false;
      symmc = 0;
      lookahead = INT_MAX;
      targetdev = -1;
      order = o;
      possibles = get_init_matrix(dim);
      selected = Matrix<int>(0, dim);
      occupation = Map<Vector<int>, int>();
      for(Entire<Rows<Matrix<int> > >::const_iterator g = entire(rows(possibles)); !g.at_end(); g++) {
         occupation[*g] = 0;
      }
   }

   bool move_no_update(Vector<int> v){
      Vector<int> third, partner;
      if(symmetric){
         partner = tof3(symmMat*v + symmShift);
      }
      if(occupation[v] != 0){
         // cout << "Returning false." << endl;
         return false;
      }
      occupation[v] = -1;
      for(Entire<Rows<Matrix<int> > >::const_iterator s = entire(rows(selected)); !s.at_end(); s++) {
         third = get_third_on_line(*s, v);
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

   int get_min_lines(){
      int result = INT_MAX;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            result = std::min(result, pt->second);
         }
      }
      return result;
   }
   
   int get_max_lines(){
      int result = INT_MIN;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            result = std::max(result, pt->second);
         }
      }
      return result;
   }
   
   int get_total_lines(){
      int result = 0;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            result += pt->second;
         }
      }
      return result;
   }
   

   int get_deviants(int minlines){
      int result = 0;
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second != -1){
            if(pt->second != minlines){
               result += pt->second == 0 ? 0 : 1;
            }
         }
      }
      return result;
   }

   Vector<int> evaluate_board(){
      // return selected.rows();
      Vector<int> result(6);
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
      Vector<int> eval = evaluate_board();
      cout << "Current evaluation:" << endl;
      cout << "Selected(0): " << eval[0] << " ";
      cout << "Possibles(5): " << eval[5] << " ";
      cout << "Deviants(2): " << eval[2] << " ";
      cout << "Minlines(1): " << eval[1] << " ";
      cout << "Maxlines(3): " << eval[3] << " ";
      cout << "TotalLines(4): " << eval[4] << " ";
      cout << endl;
   }


   Vector<int> monte_carlo_try_no_undo() {
      bool r = finish_monte_carlo();
      if(!r){ cout << "Something went wrong while playing!" << endl; }
      return evaluate_board();
   }

   Vector<int> monte_carlo_try(){
      int presize = selected.rows();
      Vector<int> result = monte_carlo_try_no_undo();
      while(selected.rows() > presize){
         undo_last_no_update();
      }
      update_possibles();
      return result;
   }

   void set_lookahead(int n){
      lookahead = n;
   }

   bool finish_monte_carlo(){
      int i = 0;
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
      possibles = Matrix<int>(0, dim);
      for(Entire<Map<Vector<int>, int> >::const_iterator pt = entire(occupation); !pt.at_end(); ++pt){
         if(pt->second == 0){
            possibles /= pt->first;
         }
      }
   }


   Vector<int> monte_carlo_evaluate(Vector<int> pt, int tries){
      int i, presize=selected.rows();
      move(pt);
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
      while(selected.rows() > presize){
         undo_last();
      }
      return result;
   }

   Vector<int> find_best_move_monte_carlo(int tries){
      Matrix<int> pts(possibles);
      // if(pts.rows() == 0){
      //    cout << endl << "I went here, even though there are no possibles." << endl;
      // }
      Matrix<int> currentBest(0,dim);
      currentBest /= possibles[0];
      Vector<int> currentBestVal = monte_carlo_evaluate(possibles[0], tries), test;
      int i=0, total = pts.rows();
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(pts)); !pt.at_end(); pt++) {
         cout << "Point no. " << i << " of " << total;
         test = monte_carlo_evaluate(*pt, tries);
         cout << " gives " << test[0] << " - ";
         // if(occupation[*pt] != 0){
         //    cout << endl << "Something went wrong! " << occupation[*pt] << endl;
         // }
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
      cout << "Winning value is: " << currentBestVal << " : " << currentBest.rows() << endl;
      i = currentBest.rows() == 1 ? 0 : rand() % currentBest.rows();
      cout << "Selecting index " << i << " : " << currentBest[i] << endl;
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
   int dim;
   Map<Vector<int>, int> pointIndices;
   Matrix<int> points, lines;

public:
   ReverseSet(int d){
      dim = d;
      points = get_init_matrix(d);
      // cout << "Points ok." << endl;
      lines(0,3);
      Vector<int> third(3), line(3);
      int i = 0, j, thirdIndex;
      int n = points.rows();
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(points)); !pt.at_end(); pt++) {
         pointIndices[*pt] = i;
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
               lines = lines/Vector<int>(line);
            }
         }
      }
   }

   ReverseSet(Matrix<int> pts, Map<Vector<int>, int> pI, Matrix<int> l){
      points = pts;
      pointIndices = pI;
      lines = l;
   }

   Matrix<int> get_points(){
      return points;
   }

   Matrix<int> get_lines(){
      return lines;
   }

   Map<Vector<int>, int> get_pt_index_map(){
      return pointIndices;
   }

   void remove_point(const Vector<int> pt){
      int pointIndex = pointIndices[pt];
      points = points.minor(~scalar2set(find_index_of_pt(pt)), All);
      Matrix<int> keepLines(0,3);
      for(Entire<Rows<Matrix<int> > >::const_iterator line = entire(rows(lines)); !line.at_end(); line++) {
         if(!line_contains(*line, pointIndex)){
            keepLines /= *line;
         }
      }
      lines = keepLines;
   }

   bool line_contains(const Vector<int> line, int index){
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

   int find_index_of_pt(const Vector<int> p){
      int i = 0;
      for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(points)); !pt.at_end(); pt++) {
         if(p == *pt){
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


int scalp(Vector<int> a, Vector<int> b){
   int result = 0, i;
   for(i=0; i<a.dim(); i++){
      result += a[i] * b[i];
   }
   return result % 3;
}

Matrix<int> points_on_hyperplane(Matrix<int> given, Vector<int> normal, int val){
   int dim = given.cols();
   Matrix<int> result(0,dim);
   for(Entire<Rows<Matrix<int> > >::const_iterator pt = entire(rows(given)); !pt.at_end(); pt++) {
      int test = scalp(*pt, normal);
      if(test == val){
         result /= *pt;
      }
   }
   return result;
}


perl::ListReturn get_init_reverse_set(int d){
   ReverseSet rs(d);
   perl::ListReturn result;
   result << rs.get_points();
   result << rs.get_pt_index_map();
   result << rs.get_lines();
   return result;
}

perl::ListReturn remove_pt_from_reverse_set(Matrix<int> pts, Map<Vector<int>, int> pI, Matrix<int> lines, Vector<int> pt){
   ReverseSet rs(pts, pI, lines);
   rs.remove_point(pt);
   perl::ListReturn result;
   result << rs.get_points();
   result << rs.get_pt_index_map();
   result << rs.get_lines();
   return result;
}

perl::ListReturn find_worst_hyperplane(Matrix<int> given){
   int dim = given.cols(), i, min = given.rows() + 1, val = 0;
   Matrix<int> hyperplanes = get_init_matrix(dim), A;
   hyperplanes = hyperplanes.minor(~scalar2set(0),All);
   Vector<int> hp;
   perl::ListReturn result;
   for(Entire<Rows<Matrix<int> > >::const_iterator n = entire(rows(hyperplanes)); !n.at_end(); n++) {
      for(i=0; i<3; i++){
         A = points_on_hyperplane(given, *n, i);
         cout << "Vector: " << *n << " val: " << i << " gives " << A.rows() << endl;
         if(A.rows() < min){
            min = A.rows();
            hp = *n;
            val = i;
         }
      }
   }
   result << hp << val << points_on_hyperplane(given, hp, val);
   return result;
}

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
   
Matrix<int> run_monte_carlo_try(Matrix<int> selected){
   SetGame set = SetGame(selected, unit_matrix<int>(6));
   // set.print();
   Vector<int> result = set.monte_carlo_try_no_undo();
   if(result[0] == 15){
      cout << "Blip!" << endl;
   }
   // set.print();
   return set.get_selected();
}

Matrix<int> run_monte_carlo_try_symmetric(Matrix<int> selected, Matrix<int> sym, Vector<int> shift){
   SetGame set = SetGame(selected, unit_matrix<int>(6));
   set.set_symmetry(sym, shift);
   // set.print();
   Vector<int> result = set.monte_carlo_try_no_undo();
   if(result[0] == 15){
      cout << "Blip!" << endl;
   }
   // set.print();
   return set.get_selected();
}

Matrix<int> finish_given_set(Matrix<int> selected, int tries, perl::OptionSet options){
   int la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<int> order = options["order"];
   SetGame set = SetGame(selected, order);
   set.set_lookahead(la);
   set.set_target_deviation(t);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<int> finish_given_set_symmetric(Matrix<int> selected, int tries, Matrix<int> A, Vector<int> b, perl::OptionSet options){
   int la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<int> order = options["order"];
   SetGame set = SetGame(selected, order);
   set.set_symmetry(A,b);
   set.set_lookahead(la);
   set.set_target_deviation(t);
   set.finish_game_monte_carlo(tries, add);
   return set.get_selected();
}

Matrix<int> finish_set(int dim, int tries, perl::OptionSet options){
   int i, la = options["lookahead"], t = options["targetdev"];
   bool add = options["add"];
   Matrix<int> order = options["order"];
   SetGame set = SetGame(dim, order);
   set.set_lookahead(la);
   set.set_target_deviation(t);
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

Vector<int> tof3(Vector<int> a){
   int i;
   Vector<int> result(a);
   for(i=0; i<result.dim(); i++){
      result[i] = (result[i]+300) % 3;
   }
   return result;
}

Matrix<int> finish_inverse_set(int dim){
   InverseSetGame IS(dim);
   return IS.finish_randomly();
}


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
