///////////////////////////////////////////////////////////////////////////////
// LRUCache.hh
// This is the LRUCache from the second answer of 
// http://stackoverflow.com/questions/2504178/lru-cache-design
// There are only minor modifications.
//
// TODO:
// 1. Improve hash function for TriangNode
///////////////////////////////////////////////////////////////////////////////

#ifndef LRUCACHE_HH
#define LRUCACHE_HH

#include <unordered_map>

template <class KEY_T, class VAL_T, class hash_f=std::hash<KEY_T> > class LRUCache{
   private:
      std::list< std::pair<KEY_T,VAL_T> > item_list;
      std::unordered_map<KEY_T, decltype(item_list.begin()), hash_f > item_map;
      size_t cache_size;
   private:
      void clean(void){
         while(item_map.size()>cache_size){
            auto last_it = item_list.end(); last_it --;
            item_map.erase(last_it->first);
            item_list.pop_back();
         }
      };
   public:
      LRUCache(int cache_size_):cache_size(cache_size_){
         ;
      };

      const std::list< std::pair<KEY_T,VAL_T> >& get_all_items(){
         return item_list;
      }

      void put(const KEY_T &key, const VAL_T &val){
         auto it = item_map.find(key);
         if(it != item_map.end()){
            item_list.erase(it->second);
            item_map.erase(it);
         }
         item_list.push_front(std::make_pair(key,val));
         item_map.insert(std::make_pair(key, item_list.begin()));
         clean();
      };
      bool exist(const KEY_T &key){
         return (item_map.count(key)>0);
      };
      const VAL_T& get(const KEY_T &key){
         assert(exist(key));
         auto it = item_map.find(key);
         item_list.splice(item_list.begin(), item_list, it->second);
         return it->second->second;
      };

};


#endif
