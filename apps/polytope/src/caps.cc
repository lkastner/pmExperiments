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
#include "polymake/group/permlib.h"
#include "polymake/polytope/LRUCache.h"



namespace polymake {
namespace polytope {
namespace {


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunused-private-field"
#pragma clang diagnostic ignored "-Wreorder"
#endif

typedef polymake::group::PermlibGroup PLG;

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

// Copypasta from 
// polymake-source/apps/group/include/named_groups.h 
Array<Array<int>>
symmetric_group_gens(int n) {
   Array<Array<int>> sgs(n-1);
   for (int i = 0; i < n-1; ++i) {
      Array<int> gen(n);
      for (int j = 0; j < n; ++j)
         gen[j] = j;
      std::swap(gen[i], gen[i+1]);
      sgs[i] = gen;
   }
   return sgs;
}

// Convert vector indicating a permutation to corresponding matrix.
Matrix<int> vector2matrix(const Vector<int> perm){
   Matrix<int> A(perm.dim(), perm.dim());
   for(int i=0; i<perm.dim(); i++){
      A(i, perm[i]) = 1;
   }
   return A;
}


class MappingData {
   private:
      int dim;
      Map<Vector<int>, int> vector2number;
      Map<int, Vector<int>> number2vector;
      Set<Vector<int>> vectors, nonZeroVectors;
      PLG plg;
      
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

      
      void init_permlib(){
         Array<Array<int>> gens(symmetric_group_gens(dim));
         plg = PLG(gens);
      }
      
   public:

      MappingData(int n): dim(n){
         init_conversions();
         init_permlib();
      }
       
       std::list<Matrix<int>> permlib_find_stabilizer(const Vector<int>& signature) const{
         cout << "Signature: " << signature << endl;
         PLG stab = plg.vector_stabilizer(signature);
         std::vector<pm::Array<int>> gp = polymake::group::all_group_elements_impl(stab);
         cout << "Permlib did its job." << endl;
         std::list<Matrix<int>> result;
         for(const auto& g: gp){
            result.push_back(vector2matrix(Vector<int>(g)));
         }
         return result;
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

      int get_dim() const{
         return dim;
      }

      int get_npoints() const {
         return vectors.size();
      }
   
};

class Cap {
   private:
      const MappingData& md;
      Set<int> occupied;
      Vector<int> linesThroughPoints, pointsOnHyperplanes;

   public:

      Cap(const MappingData& md_in): 
         md(md_in),
            linesThroughPoints(md_in.get_npoints()), pointsOnHyperplanes(3*md_in.get_npoints()){
         }

      Cap(Set<int> o,
         Vector<int> ltp,
         Vector<int> poh,
         const MappingData& md_in): occupied(o), linesThroughPoints(ltp), pointsOnHyperplanes(poh), md(md_in){
         }

      Cap mutate(const Matrix<int>& A, const Vector<int>& b){
         Cap result(md);
         Vector<int> s;
         for(const auto& e: occupied){
            s = md.get_number2vector()[e];
            s = (A*s) + b;
            vectorMod3(s);
            result.select(s);
         }
         return result;
      }

      Vector<int> hyperplane_indicator(const Vector<int>& h) const{
         Vector<int> result(3);
         for(const auto& o: occupied){
            int eval = h * md.get_number2vector()[o];
            eval %= 3;
            eval += 3;
            eval %= 3;
            result[eval]++;
         }
         // std::sort(result.begin(), result.end());
         // std::reverse(result.begin(), result.end());
         return result;
      }

      void select(const Vector<int>& v){
         for(const auto& j: occupied){
            const Vector<int>& third = get_third_on_line(v, md.get_number2vector()[j]);
            linesThroughPoints[md.get_vector2number()[third]]++;
         }
         // for(int i=0; i<md.get_npoints(); i++){
         //    int val = (md.get_number2vector()[i]) * v;
         //    val %= 3;
         //    val += 3;
         //    val %= 3;
         //    pointsOnHyperplanes[3*i+val]++;
         // }
         occupied += md.get_vector2number()[v];
      }

      void select(const int n){
         select(md.get_number2vector()[n]);
      }
      
      void unselect(const Vector<int>& v){
         occupied -= md.get_vector2number()[v];
         for(const auto& j: occupied){
            const Vector<int>& third = get_third_on_line(v, md.get_number2vector()[j]);
            linesThroughPoints[md.get_vector2number()[third]]--;
         }
         // for(int i=0; i<md.get_npoints(); i++){
         //    int val = (md.get_number2vector()[i]) * v;
         //    val %= 3;
         //    val += 3;
         //    val %= 3;
         //    pointsOnHyperplanes[3*i+val]--;
         // }
      }

      void print_hyperplane_arrangement() const{
         for(int i=0; i<md.get_npoints(); i++){
            const Vector<int>& hyp = (md.get_number2vector()[i]);
            cout << hyp << ": ";
            Vector<int> result(hyperplane_indicator(hyp));
            cout << result << endl;
         }
      }
      
      void unselect(const int n){
         unselect(md.get_number2vector()[n]);
      }

      bool checkBound(const int& b){
         // Something is wrong here.
         for(const auto& e: pointsOnHyperplanes){
            if(e>b){
               print();
               cout << e << endl;
               cout << "Bound violated." << endl;
               return false;
            }
         }
         return true;
      }

      void print() const{
         for(const auto& e: occupied){
            cout << e << " ";
         }
         cout << endl;
         for(const auto& e: occupied){
            cout << (md.get_number2vector()[e]) << ", ";
         }
         cout << endl;
         // for(const auto& e: linesThroughPoints){
         //    cout << e << " ";
         // }
         // cout << endl;
         cout << signature() << endl;
         cout << endl;
      }

      std::map<int, Set<int>> filterLTP() const{
         std::map<int, Set<int>> result;
         for(int i=0; i<linesThroughPoints.dim(); i++){
            result[-linesThroughPoints[i]] += i;
         }
         return result;
      }

      void batch_select(const Set<int>& o){
         occupied.clear();
         linesThroughPoints = Vector<int>(md.get_npoints());
         for(const auto& no: o){
            select(no);
         }
      }

      Vector<int> get_ith_filter(const int i) const{
         Vector<int> result(3);
         for(const auto& o : occupied){
            result[md.get_number2vector()[o][i]]++;
         }
         return result;
      }

      bool isValid() const{
         Vector<int> v;
         for(int i=0; i<md.get_dim(); i++){
            v = get_ith_filter(i);
            if(v[0] < v[1]) return false;
            if(v[1] < v[2]) return false;
         }
         return true;
      }

      Set<int> get_free_points() const{
         Set<int> result;
         for(int i=0; i<linesThroughPoints.dim(); i++){
            if(linesThroughPoints[i] == 0){
               result += i;
            }
         }
         result -= occupied;
         return result;
      }
      
      Array<Vector<int>> signature() const{
         Array<Vector<int>> result(md.get_dim());
         for(int i=0; i<md.get_dim(); i++){
            result[i] = get_ith_filter(i);
         }
         return result;
      }

      Vector<int> signature_vector() const{
         Array<Vector<int>> sig(signature());
         Vector<int> result(sig.size());
         int counter=0;
         for(int i=1; i<sig.size(); i++){
            if(sig[i] != sig[i-1]){
               counter++;
            }
            result[i] = counter;
         }
         return result;
      }

      Vector<int> get_permutation() const{
         Vector<int> start(range(0,md.get_dim()-1));
         Array<Vector<int>> F(signature());
         // cout << "Before: " << F << endl;
         std::sort(start.begin(), start.end(),[&](const int& a, const int& b) {
            return lex_compare(F[a], F[b]) == -1;
         });
         // cout << "After: " << start << endl;
         return start;
      }

      Set<int> apply_matrix_to_pts(const Matrix<int> A) const{
         Set<int> newp;
         for(const auto& o: occupied){
            newp += md.get_vector2number()[A*(md.get_number2vector()[o])];
         }
         return newp;
      }

      void apply_matrix_transform(const Matrix<int> A){
         batch_select(apply_matrix_to_pts(A));
      }

      void canonicalize() {
         Vector<int> perm = get_permutation();
         Matrix<int> A(vector2matrix(perm));
         // cout << "Transform: " << endl;
         // cout << A << endl;
         apply_matrix_transform(A);
         // std::list<Matrix<int>> stab(md.permlib_find_stabilizer(signature_vector()));
         // cout << "Stab ok." << endl;
         // Map<Set<int>, Matrix<int>> sorted;
         // for(const auto& B : stab){
         //    cout << "B: " << B << endl;
         //    Set<int> key(apply_matrix_to_pts(B));
         //    sorted[key] = B;
         // }
      }

      Cap predecessor() const{
         Cap result(md);
         result.batch_select(get_content());
         result.canonicalize();
         // cout << "Pred: " << result.get_content() << endl;
         // cout << "Last: " << result.get_content().back() << endl;
         int back = result.get_content().back();
         result.unselect(back);
         result.canonicalize();
         return result;
      }

      const Set<int>& get_content() const{
         return occupied;
      }

      // const MappingData& get_md() const{
      //    return md;
      // }

      // Cap(const Cap &c):md(c.get_md()){
      //    batch_select(c.get_content());
      // }

      Cap & operator= (const Cap & other){
         if (this != &other){ // protect against invalid self-assignment
            this->batch_select(other.get_content());
         }
         return *this;
      }

      std::vector<Cap> get_children() const{
         std::vector<Cap> result;
         Set<Array<Vector<int>>> signatures;
         int last = occupied.back();
         const Set<int> free = get_free_points();
         for(const auto& f:free){
            if(f > last){
               // cout << "f: " << f << " - " << (md.get_number2vector()[f]) << endl;
               Cap candidate(md);
               candidate.batch_select(get_content());
               candidate.select(f);
               candidate.canonicalize();
               if(signatures.contains(candidate.signature())){
                  // cout << "Hello." << endl;
               } else {
                  signatures += candidate.signature();
                  result.push_back(candidate);
                  // if(candidate.isValid()){
                  //    result.push_back(candidate);
                  // } else {
                  //    // cout << "Invalid candidate." << endl;
                  // }
               }
            }
         }
         return result;
      }

      int size() const{
         return occupied.size();
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
      const Set<Vector<int>>& vectors, nonZeroVectors;
      Map<int, Set<int>> supports;
      

   public:
      Optimizer(int n, const MappingData& md) : dim(n),
         vector2number(md.get_vector2number()),
         number2vector(md.get_number2vector()),
         vectors(md.get_vectors()),
         nonZeroVectors(md.get_nonZeroVectors())
         {
         Matrix<int> m({{2,1,1},{0,2,0},{0,0,1}});
         f3inverse(m);
      }
      
};

class ReverseSearch {
   private:
      int dim, max;
      const MappingData md;
      Cap currentMax;
      LRUCache<Set<int>, std::vector<Cap>, pm::hash_func<Set<int>>> neighborCache;

      
      const std::vector<Cap>& get_children(const Cap& c){
         if(!neighborCache.exist(c.get_content())){
            neighborCache.put(c.get_content(), c.get_children());
         }
         return neighborCache.get(c.get_content());
      }
   
      int nChildren(const Cap& c){
         // TODO: Cache this.
         // cout << "nchildren. " << (c.get_content()) << endl;
         return get_children(c).size();
      }

      Cap predecessor(const Cap& c){
         // TODO: Cache this.
         // cout << "pred." << endl;
         return c.predecessor();
      }

      Cap jthChild(const Cap& c, int j){
         // TODO: Cache this.
         // cout << "jthchild. j: " << j << endl;
         const std::vector<Cap>& neighbors(get_children(c));
         // cout << (neighbors.size()) << endl;
         return neighbors[j];
      }

      void numbered_predecessor(Cap& c, int& j){
         // Set c to be predecessor and j to be the number of (old) c as child
         // of (new) c.
         // cout << "Numbered pred." << endl;
         Cap pred = predecessor(c);
         const std::vector<Cap>& neighbors(get_children(pred));
         for(int i=0; i<neighbors.size(); i++){
            if(neighbors[i] == c){
               c = pred;
               j = i;
               return;
            }
         }
      }


   public:
      ReverseSearch(int d): md(d), currentMax(md), neighborCache(5000){
         max = 0;
      }
      
      void generic_reverse_search(Cap start){
         Cap v = start;
         int j=-1, depth=0;
         while(!((depth==0) && (j==nChildren(v)-1))){
            while(j<nChildren(v)-1){
               j++;
               Cap Avj = jthChild(v,j);
               if(predecessor(Avj) == v){
                  if(v.get_free_points().size() + v.size() > max){
                     v = Avj;
                     if(v.size() % 20 == 0){
                        cout << "Descending. " << v.size() << endl;
                     }
                     j = -1;
                     depth++;
                  }
               }
            }
            if(v.size() > max){
               max = v.size();
               cout << "Max: " << max << endl;
               currentMax = v;
               v.print();
            }
            if(depth > 0){
               numbered_predecessor(v, j);
               depth--;
            }
         }
         cout << "Max was: " << max << endl;
         currentMax.print();
      }
};

void do_stuff(int i){
   MappingData m(i);
   Optimizer o(i, m);
}

void do_stuff_matrix(Matrix<int> m){
   int i = m.cols();
   MappingData md(i);
   cout << "MD done." << endl;
   Cap C(md);
   for(const auto& row: rows(m)){
      C.select(row);
   }
   cout << "Printing" << endl;
   C.print();
   cout << "Done Printing" << endl;
   for(int k =0; k<i; k++){
      cout << k << ": " << C.get_ith_filter(k) << endl;
   }
   cout << "This worked." << endl;
   cout << C.get_permutation() << endl;
   cout << "Permutation done." << endl;
   C.print();
   cout << "Canonicalizing." << endl;
   C.canonicalize();
   C.print();
   cout << "Neighbors: " << (C.get_children().size())<<endl;
   // for(const auto& n: C.get_children()){
   //    n.print();
   // }
   C.predecessor().print();
   cout << "Starting RS." << endl;
   ReverseSearch RS(i);
   RS.generic_reverse_search(C);
}

void analyze_cap(Matrix<int> m){
   int i = m.cols();
   MappingData md(i);
   cout << "MD done." << endl;
   Cap C(md);
   for(const auto& row: rows(m)){
      C.select(row);
   }
   cout << "Printing" << endl;
   C.print();
   C.print_hyperplane_arrangement();
   
}

perl::ListReturn mutate_cap(Matrix<int> m){
   int i = m.cols();
   MappingData md(i);
   cout << "MD done." << endl;
   Cap C(md);
   for(const auto& row: rows(m)){
      C.select(row);
   }
   cout << "Printing" << endl;
   C.print();
   perl::ListReturn result;
   result << Vector<int>(C.get_content());
   
   for(const auto& v: md.get_vectors()){
      Cap C1 = C.mutate(unit_matrix<int>(i), v);
      Set<int> check = (C1.get_content() * C.get_content());
      if(check.size() == 0){
         result << Vector<int>(v);
         cout << v << endl;
      }
   }

   return result;
   
}

bool disjoint_cap(Matrix<int> m, const Array<Vector<int>> shifts){
   int i = m.cols();
   MappingData md(i);
   Cap C(md);
   for(const auto& row: rows(m)){
      C.select(row);
   }
   perl::ListReturn result;
   result << Vector<int>(C.get_content());
   
   std::list<Set<int>> caps;
   caps.push_back(C.get_content());
   for(const auto& v: shifts){
      Cap C1 = C.mutate(unit_matrix<int>(i), v);
      Cap C2 = C.mutate(unit_matrix<int>(i), 2*v);
      for(const auto& pc: caps){
         Set<int> check1 = (C1.get_content() * pc);
         Set<int> check2 = (C2.get_content() * pc);
         if(check1.size() != 0){
            cout << "Fail: " << v << endl;
            return false;
         }
         if(check2.size() != 0){
            cout << "Fail: " << (2*v) << endl;
            return false;
         }
      }
      Set<int> check12 = (C1.get_content() * C2.get_content());
      if(check12.size() == 0){
         caps.push_back(C1.get_content());
         caps.push_back(C2.get_content());
      } else {
         cout << "Fail12: " << v << endl;
         return false;
      }
   }
   return true;
   
}


#if defined(__clang__)
#pragma clang diagnostic pop
#endif
} // namespace

Function4perl(&do_stuff, "do_stuff");

Function4perl(&do_stuff_matrix, "do_stuff_matrix");

Function4perl(&analyze_cap, "analyze_cap");

Function4perl(&mutate_cap, "mutate_cap");

Function4perl(&disjoint_cap, "disjoint_cap");

} // namespace polymake
} // namespace polytope
