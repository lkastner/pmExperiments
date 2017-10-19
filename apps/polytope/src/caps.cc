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

Matrix<int> f3inverse(Matrix<int> in){
   Matrix<Rational> tmp(in);
   tmp = inv(tmp);
   tmp = det(in) * tmp;
   if((det(in) + 1) % 3 == 0){
      tmp = (-1) * tmp;
   }
   Matrix<int> result(in);
   for(int i=0; i<result.rows(); i++){
      for(int j=0; j<result.cols(); j++){
         int converted = pm::convert_to<int>(tmp(i,j));
         converted %= 3;
         converted += 3;
         result(i,j) = converted % 3;
      }
   }
   int testdet = det(result * in);
   if((testdet - 1) % 3 != 0){
      cout << "Something went wrong" << endl;
      cout << in << endl;
      cout << result << endl;
   }
   return result;
}

void vectorMod3(Vector<int>& b){
   for(int i=0; i<b.dim(); i++){
      b[i] %= 3;
      b[i] += 3;
      b[i] %= 3;
   }
}

Vector<int> get_third_on_line(Vector<int> a, Vector<int> b){
   Vector<int> result = -a-b;
   vectorMod3(result);
   return result;
}


class MappingData {
   private:
      int dim;
      Map<Vector<int>, int> vector2number;
      Map<int, Vector<int>> number2vector;
      Set<Vector<int>> vectors, nonZeroVectors;
      
      void get_vectors(Vector<int> current, int pos){
         if(pos == dim){
            vectors += Vector<int>(current);
            if(current != zero_vector<int>(dim)){
               nonZeroVectors += Vector<int>(current);
            }
         } else {
            for(int i=0; i<3; i++){
               current[pos] = i;
               get_vectors(current, pos+1);
            }
         }
      }

      void init_conversions(){
         Vector<int> start(dim);
         get_vectors(start, 0);
         Set<Vector<int>> tmp;
         Matrix<int> todp(dim,dim), fromdp(dim, dim);
         todp(0, 0) = 1;
         fromdp(0, 0) = 1;
         for(int i=1; i<dim; i++){
            todp(0,i) = 1;
            fromdp(0,i) = -1;
            todp(i, dim-i) = 1;
            fromdp(i, dim-i) = 1;
         }
         // cout << todp << endl;
         // cout << fromdp << endl;
         // cout << (todp * fromdp) << endl;
         for(auto v: vectors){
            tmp += todp * v;
         }
         int i=0;
         for(auto v: tmp){
            Vector<int> key(fromdp * v);
            vector2number[key] = i;
            number2vector[i] = key;
            // cout << i << ": " << key << endl;
            i++;
         }
      }

   public:
      MappingData(int n): dim(n){
         init_conversions();
      }

      const Map<Vector<int>, int>& get_vector2number() const{
         return vector2number;
      }  

      const Map<int, Vector<int>>& get_number2vector() const{
         return number2vector;
      }

      const Set<Vector<int>>& get_vectors() const{
         return vectors;
      }
      
      const Set<Vector<int>>& get_nonZeroVectors() const{
         return nonZeroVectors;
      }
   
};

class Cap {
   private:
      Set<int> occupied;
      Vector<int> linesThroughPoints;
      const Map<Vector<int>, int>& vector2number;
      const Map<int, Vector<int>>& number2vector;

   public:
      Cap(const Map<Vector<int>, int>& v2n,
      const Map<int, Vector<int>>& n2v):
         vector2number(v2n), number2vector(n2v),
            linesThroughPoints(v2n.size()){
         }

      Cap(const MappingData& md): 
         vector2number(md.get_vector2number()), number2vector(md.get_number2vector()),
            linesThroughPoints(md.get_vector2number().size()){
         }

      Cap(Set<int> o,
         Vector<int> ltp,
         const MappingData& md): occupied(o), linesThroughPoints(ltp), 
         vector2number(md.get_vector2number()), number2vector(md.get_number2vector()){
         }

      Cap mutate(const Matrix<int>& A, const Vector<int>& b){
         Cap result(vector2number, number2vector);
         Vector<int> s;
         for(const auto& e: occupied){
            s = number2vector[e];
            s = (A*s) + b;
            vectorMod3(s);
            result.select(s);
         }
         return result;
      }

      Vector<int> hyperplane_indicator(const Vector<int>& h) const{
         Vector<int> result(3);
         for(const auto& o: occupied){
            int eval = h * number2vector[o];
            eval %= 3;
            eval += 3;
            eval %= 3;
            result[eval]++;
         }
         // std::sort(result.begin(), result.end());
         // std::reverse(result.begin(), result.end());
         return result;
      }

      Cap get_subcap(const Vector<int>& h, const int& val) const {
         return Cap(vector2number, number2vector);
      }

      void select(const Vector<int>& v){
         for(const auto& j: occupied){
            const Vector<int>& third = get_third_on_line(v, number2vector[j]);
            linesThroughPoints[vector2number[third]]++;
         }
         occupied += vector2number[v];
      }
      
      void unselect(const Vector<int>& v){
         occupied -= vector2number[v];
         for(const auto& j: occupied){
            const Vector<int>& third = get_third_on_line(v, number2vector[j]);
            linesThroughPoints[vector2number[third]]--;
         }
      }

      void print(){
         for(const auto& e: occupied){
            cout << e << " ";
         }
         cout << endl;
         for(const auto& e: linesThroughPoints){
            cout << e << " ";
         }
         cout << endl;
      }

      std::map<int, Set<int>> filterLTP() const{
         std::map<int, Set<int>> result;
         for(int i=0; i<linesThroughPoints.dim(); i++){
            result[-linesThroughPoints[i]] += i;
         }
         return result;
      }

      const Set<int>& get_content() const{
         return occupied;
      }
      
};
      
bool operator==(const Cap& c1, const Cap& c2){
   return c1.get_content() == c2.get_content();
}

class Optimizer {
   private:
      int dim, multCount;
      const Map<Vector<int>, int>& vector2number;
      const Map<int, Vector<int>>& number2vector;
      Map<int, Map<int, std::pair<Matrix<int>, Vector<int>>>> switchTable;
      const Set<Vector<int>>& vectors, nonZeroVectors;
      Map<int, Set<int>> supports;


      void init_switchTable(){
         Matrix<int> I(unit_matrix<int>(dim));
         for(auto v: vectors){
            Vector<int> b(-v);
            vectorMod3(b);
            switchTable[0][vector2number[v]] = std::pair<Matrix<int>, Vector<int>>(I, b);
            supports[0] += vector2number[v];
         }
         Vector<int> b(zero_vector<int>(dim));
         for(int i=1; i<=dim; i++){
            for(auto v: vectors){
               Matrix<int> tmp(unit_matrix<int>(i-1) / zero_matrix<int>(dim-i+1, i-1));
               tmp = tmp | v;
               tmp = tmp | (zero_matrix<int>(i-1, dim-i+1) / unit_matrix<int>(dim-i+1));
               // cout << tmp << endl;
               for(int j=i; j<dim+1; j++){
                  int test = det(tmp.minor(All, ~scalar2set(j)));
                  if(test != 0){
                     Matrix<int> A(tmp.minor(All, ~scalar2set(j)));
                     A = f3inverse(A);
                     switchTable[i][vector2number[v]] = std::pair<Matrix<int>, Vector<int>>(A, b);
                     supports[i] += vector2number[v];
                     break;
                  }
               }
            }
         }
      }

      void print_switchTable(){
         Integer totalsize = 1;
         for(auto& row: switchTable){
            int i = row.first;
            cout << "Support: {";
            for(auto& e: supports[i]){
               cout << e << " ";
            }
            cout << "}" << endl;
            totalsize *= supports[i].size();
            for(auto& entry: row.second){
               int j = entry.first;
               cout << "i: " << i << " - " << number2vector[i] << endl;
               cout << "j: " << j << " - " << number2vector[j] << endl;
               Matrix<int>& A(entry.second.first);
               Vector<int>& b(entry.second.second), tmp(A*number2vector[j]);
               tmp += b;
               vectorMod3(tmp);
               cout << "b: " << b<< endl;
               cout << "result: " << tmp << endl;
               cout << A << endl;
            }
         }
         cout << "Total group size: " << totalsize << endl;
      }

      const std::list<std::pair<Matrix<int>, Vector<int>>> get_switches(Set<int> support, int fixed) const{
         std::list<std::pair<Matrix<int>, Vector<int>>> result;
         Set<int> intersection = support * supports[fixed];
         for(auto& e: intersection){
            result.push_back(switchTable[fixed][e]);
         }
         return result;
      }

      Cap optimize(const Cap& in, int fixed){
         if(fixed == dim+1){
            return Cap(in);
         }
         Cap result(in);
         std::map<int, Set<int>> LTP = in.filterLTP();
         std::map<int, Set<int>>::const_iterator desired = LTP.begin();
         std::list<std::pair<Matrix<int>, Vector<int>>> switches = get_switches(desired->second, fixed);
         while(switches.size() == 0){
            ++desired;
            switches = get_switches(desired->second, fixed);
         }
         for(const auto& s: switches){
            optimize(result.mutate(s.first, s.second), fixed+1);
            multCount++;
            if(multCount % 1000 == 0){
               cout << "Took " << multCount << " multiplications so far." << endl;
            }
         }
         return result;
      }
      

   public:
      Optimizer(int n, const MappingData& md) : dim(n),
         vector2number(md.get_vector2number()),
         number2vector(md.get_number2vector()),
         vectors(md.get_vectors()),
         nonZeroVectors(md.get_nonZeroVectors())
         {
         init_switchTable();
         Matrix<int> m({{2,1,1},{0,2,0},{0,0,1}});
         f3inverse(m);
         print_switchTable();
      }
      
      Cap optimize(const Cap& in){
         multCount = 0;
         Cap result = optimize(in, 0);
         cout << "Took " << multCount << " multiplications." << endl;
         return result;
      }

      void optimize2(const Cap& in){
         Map<Vector<int>, int> tmp;
         for(const auto& v: nonZeroVectors){
            tmp[in.hyperplane_indicator(v)]++;
         }
         for(const auto& p: tmp){
            cout << p.first << ": " << p.second << endl;
         }
      }
};

class ReverseSearch {
   private:
      int dim;
      const MappingData md;
      const Optimizer opt;

      int nChildren(Cap c){
         // Get number of neighbors of c
         return 10;
      }

      Cap predecessor(const Cap& c){
         // Get predecessor of c
         return Cap(c);
      }

      Cap jthChild(const Cap& c, int j){
         // Get the jth neighbor of c
         return Cap(c);
      }

      void numbered_predecessor(Cap& c, int& j){
         // Set c to be predecessor and j to be the number of (old) c as child
         // of (new) c.
      }

      void generic_reverse_search(Cap start){
         Cap v = start;
         int j=-1, depth=0;
         while(!((depth==0) && (j==nChildren(v)))){
            while(j<nChildren(v)){
               j++;
               Cap Avj = jthChild(v,j);
               if(predecessor(Avj) == v){
                  //v = Avj;
                  j = -1;
                  depth++;
               }
            }
            if(depth > 0){
               numbered_predecessor(v, j);
               depth--;
            }
         }
      }

   public:
      ReverseSearch(int d): md(d), opt(d, md){
      }
};

void do_stuff(int i){
   MappingData m(i);
   Optimizer o(i, m);
}

void do_stuff_matrix(Matrix<int> m){
   int i = m.cols();
   MappingData md(i);
   Optimizer o(i, md);
   cout << "Optimizer done." << endl;
   Cap C(md);
   for(const auto& row: rows(m)){
      C.select(row);
   }
   C.print();
   cout << "Will optimize now." << endl;
   o.optimize2(C);
}


#if defined(__clang__)
#pragma clang diagnostic pop
#endif
} // namespace

Function4perl(&do_stuff, "do_stuff");

Function4perl(&do_stuff_matrix, "do_stuff_matrix");

} // namespace polymake
} // namespace polytope
