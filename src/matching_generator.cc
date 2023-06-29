#include "matching_generator.hh"
#include "iostream"
#include <stdarg.h>
#include <stack>
#include <unordered_set>
#include <sys/time.h>
using namespace std;

MatchingGenerator::MatchingGenerator( uint8_t stripe_x, uint16_t _node_num, uint8_t k, uint8_t m,
                                     uint _stripes_num)
    : x{stripe_x}, node_num{_node_num}, rs_k{k}, rs_m{m}, stripes_num{_stripes_num} {
  rs_n = (x-1) * rs_k + (x-1) * rs_m;
  Stripe::set_params(node_num, rs_k, rs_m);
  stripes.resize(stripes_num);
  //size_t table_size = static_cast<size_t>(stripes_num);
  //table_size = table_size * (table_size - 1) / 2;
  cost_table.assign(stripes_num / x, 0);
  //table = build_search_table(x);
  if(x == 2){
    table = {
      {0,{{1,0,0,0}}},
      {1,{{0,1,0,0}}},
      {2,{{0,0,1,0}}},
      {3,{{0,0,0,1}}}
    };
  }
  if(x == 3){
    table = {
      {0,{{2,0,0,0}}},
      {1,{{1,1,0,0}}},
      {2,{{1,0,1,0},
          {0,2,0,0}}},
      {3,{{1,0,0,1},
          {0,1,1,0}}},
      {4,{{0,1,0,1},
          {0,0,2,0}}},
      {5,{{0,0,1,1}}},
      {6,{{0,0,0,2}}}
    };
  }
  if(x == 4){
    table = {
      {0,{{3,0,0,0}}},
      {1,{{2,1,0,0}}},
      {2,{//{2,0,1,0},
          {1,2,0,0}}},
      {3,{{2,0,0,1},
          {1,1,1,0},
          {0,3,0,0}}},
      {4,{//{1,1,0,1},
          {1,0,2,0},
          {0,2,1,0}
                   }},
      {5,{{1,0,1,1},
          {0,1,2,0},
          //{0,2,0,1}
                   }},
      {6,{{1,0,0,2},
          {0,0,3,0},}},
      {7,{{0,1,0,2}}},
      {8,{{0,0,1,2}}},
      {9,{{0,0,0,3}}}
    };
  }
}
std::unordered_map<uint, set<std::vector<int>>> MatchingGenerator::build_search_table(int X){
  if(X == 2){
    std::unordered_map<uint, set<std::vector<int>>> temptable;
    for(int i = 0; i < rs_m ; i++){
      vector<int> search_vector(rs_m, 0);
      search_vector[i] = 1;
      temptable[i].insert(search_vector);
    }
    return temptable;
  }
  else{
    std::unordered_map<uint, set<std::vector<int>>> front_temptable = build_search_table(X-1);
    std::unordered_map<uint, set<std::vector<int>>> temptable;
    for(int i = 0;i <= (X-1)*(rs_m - 1); i++){
      for(int j = 0;j <= i && j < rs_m;j++){
        if(front_temptable.find(i-j) == front_temptable.end()) continue;
        if(front_temptable[i-j].size()){
          for(auto it : front_temptable[i-j]){
            vector<int> search_vector = it;
            search_vector[j]++;
            temptable[i].insert(search_vector);
          }
        }
      }
    }
    return temptable;
  }
}
void MatchingGenerator::build_map() {
  const uint8_t upside = 1 << rs_m;
  std::u16string key;
  parity_map.resize(K);
  partial_dist_map.resize(K);
  for(uint cluster = 0; cluster < K ; cluster++){
    auto map_it = parity_map[cluster].end();
    auto index_it = partial_dist_map[cluster].end();
    
    for (uint index = 0; index < parity_cluster_0[cluster].size(); ++index) {
      const std::u16string &pd = stripes[parity_cluster_0[cluster][index]].parity_dist;
      // map parity dist to its correspondent stripes' indexes
      if ((map_it = parity_map[cluster].find(pd)) == parity_map[cluster].end()) {
        parity_map[cluster][pd] = std::list<uint>{parity_cluster_0[cluster][index]};
      } else {
        map_it->second.push_back(parity_cluster_0[cluster][index]);
      }
      
      // build multiple index for parity dist
      
      for (uint8_t s = 1, bits; s < upside; ++s) {
        key.clear();
        bits = s;
        for (auto &x : pd) {
          if (bits & 1) {
            key.push_back(x);
          }
          bits >>= 1;
        }
        key.push_back(s);
        if ((index_it = partial_dist_map[cluster].find(key)) == partial_dist_map[cluster].end()) {
          partial_dist_map[cluster][key] = std::list<uint>{parity_cluster_0[cluster][index]};
        } else {
          index_it->second.push_back(parity_cluster_0[cluster][index]);
        }
      }
    }
  }
}
void MatchingGenerator::mark_matched(uint i, std::vector<uint> j, uint8_t cost, bool flag) {
  stripes[i].matching_index = j;
  stripes[i].selected = true;
  stripes[i].cur_cost = cost;
  for(uint it = 0 ;it < j.size() ; it++){
    std::vector<uint> temp_vetctor = j;
    temp_vetctor[it] = i;
    stripes[j[it]].matching_index = temp_vetctor;
    stripes[j[it]].selected = true;
    stripes[j[it]].cur_cost = cost;
  }
  if (flag) {
    if (!cost) {
      const std::u16string &pd = stripes[i].parity_dist;
      auto map_it = zero_map.find(pd);
      if (map_it != zero_map.end()) {
        map_it->second.push_front(i);
      } else {
        zero_map[pd] = std::list<uint>{i};
      }
    } else{
      bool a = true;
      for(auto it : j){
        if (stripes[it].parity_dist != stripes[it].parity_dist) a = false;
      }
      if (a) {
        none_zero_candidate.push_back(i);
      }
    }
  }
}

bool MatchingGenerator::single_search(uint8_t target_cost, uint index, uint temp_index, 
                    std::list<uint>::iterator temp_it,  std::vector<int> search_vector, 
                    std::vector<std::list <uint>> search_domain, std::vector<bool> flag, std::vector<int> temp_count,
                    uint total_conut,
                    std::vector<uint> tempvector,
                    std::queue<uint> cleanup_queue,
                    uint8_t min_cost,
                    std::vector<uint> min_cost_index) {
while(search_vector[temp_index] == 0){
  temp_index++;
}
if(temp_index >= rs_m) return false;
if(search_domain[temp_index].size() < search_vector[temp_index]) return false;
std::vector<std::list<uint>::iterator> it_vector;
temp_it = search_domain[temp_index].begin();
while(temp_it != search_domain[temp_index].end()){
  if(!stripes[*temp_it].selected){
    it_vector.emplace_back(temp_it);
    flag[*temp_it] = true;
    temp_count[temp_index]++;
    total_conut++;
    break;
  }
  else{
    temp_it++;
  }
}
if(temp_it == search_domain[temp_index].end()) return false;
while(it_vector.size()){
  tempvector.clear();
  tempvector.emplace_back(index);
  for(auto it : it_vector){
    tempvector.emplace_back(*(it));
  }
  if(total_conut == x) {
    if(get_cost(tempvector) < min_cost){
      min_cost = get_cost(tempvector);
      min_cost_index.assign(tempvector.begin() + 1 , tempvector.end());
      if(min_cost == target_cost) break;
      else{
        temp_it = it_vector.back();
        it_vector.pop_back();
        flag[*temp_it] = false;
        temp_count[temp_index]--;
        total_conut--;
      }
    }
    else{
      temp_it = it_vector.back();
      it_vector.pop_back();
      flag[*temp_it] = false;
      temp_count[temp_index]--;
      total_conut--;
    }
  }
  else{
    if(get_cost(tempvector) > target_cost){
      temp_it = it_vector.back();
      it_vector.pop_back();
      temp_count[temp_index]--;
      total_conut--; 
      flag[*temp_it] = false;
    }
  }
  if(temp_count[temp_index] < search_vector[temp_index]){
    while(++temp_it != search_domain[temp_index].end()){
      if(!flag[*temp_it] && !stripes[*temp_it].selected){
        flag[*temp_it] = true;
        it_vector.emplace_back(temp_it);
        temp_count[temp_index]++;
        total_conut++;
        break;
      }
    }
    if(temp_it != search_domain[temp_index].end());
    else{
      if(!it_vector.size())return false;
      if(temp_count[temp_index] == 0 && temp_index > 0){
        temp_index--;
        temp_it = it_vector.back();
        it_vector.pop_back();
        temp_count[temp_index]--;
        total_conut--;
        flag[*temp_it] = false;
      }
      else{
        temp_it = it_vector.back();
        it_vector.pop_back();
        temp_count[temp_index]--;
        total_conut--;
        flag[*temp_it] = false;
      }
    }
  }
  else{
    temp_index++;
    if(temp_index >= rs_m) return false;
    temp_it = search_domain[temp_index].begin();
    while(temp_it != search_domain[temp_index].end()){
      if(!flag[*temp_it] && !stripes[*temp_it].selected){
        it_vector.emplace_back(temp_it);
        flag[*temp_it] = true;
        temp_count[temp_index]++;
        total_conut++;
        break;
      }
      else{
        temp_it++;
      }
    }
  }
}

if(min_cost < stripes[index].cur_cost){
  stripes[index].cur_cost = min_cost;
  stripes[index].matching_index = min_cost_index;
  if(min_cost == target_cost){
    cleanup_queue.push(index);
    if(it_vector.size()){
      for(auto it : min_cost_index){
        cleanup_queue.push(it);
      }
    }
  }
}
return min_cost == target_cost;
}  

bool MatchingGenerator::search_for_matching(uint index, uint8_t target_cost, uint cluseter) {
  static std::queue<std::list<std::u16string>::iterator> cleanup_queue;
  static std::u16string key;
  static const uint8_t all_one = (1 << rs_m) - 1;
  const std::u16string &pd = stripes[index].parity_dist;
  std::vector<std::list <uint>> search_domain;
  search_domain.resize(rs_m);
  if(parity_map[cluseter][pd].size() > 1){
    for(auto it = parity_map[cluseter][pd].begin(); it != parity_map[cluseter][pd].end(); it++){
      if(*it != index) search_domain[0].push_back(*it);
    }
  }
  for(int i = 1 ; i < search_domain.size(); i++){
    int x, y, bits;
    int k = static_cast<int>(rs_m - i), comb = (1 << k) - 1;
    auto index_it = partial_dist_map[cluseter].end();
    // black magic to get k-subset
    while (comb < all_one) {
      key.clear();
      bits = comb;
      for (auto &x : pd) {
        if (bits & 1) {
          key.push_back(x);
        }
        bits >>= 1;
      }
      key.push_back(static_cast<char16_t>(comb));
      if ((index_it = partial_dist_map[cluseter].find(key)) != partial_dist_map[cluseter].end() &&
          index_it->second.size() > 1) {
          for(auto it = index_it->second.begin(); it != index_it->second.end(); it++){
            if(*it != index) search_domain[i].push_back(*it);
          }
      }
      // get next comb of k-subset
      x = comb & (-comb), y = x + comb;
      comb = (((comb & ~y) / x) >> 1) | y;
    }
    search_domain[i].sort();
    search_domain[i].unique();
  }
  /*
  cout <<"stripe: "<< index << "'s search domain " << endl;
  for(uint i = 0; i < search_domain.size(); i++){
    cout << i << " :";
    for(auto it : search_domain[i])
      cout << it << " ";
    cout << endl;
  }
  */
  bool flag_select = false;
  for(auto search_vector:table[target_cost]){
    //std::vector<int> search_vector = table[target_cost][i];
    std::vector<bool> flag(stripes_num, 0);
    std::vector<int> temp_count(rs_m, 0);
    uint temp_index = 0;
    uint total_conut = 1;
    std::vector<uint> tempvector = {index};
    static std::queue<uint> cleanup_queue;
    uint8_t min_cost = stripes[index].cur_cost;
    uint8_t tmp_cost;
    std:;vector<uint> min_cost_index = {UINT32_MAX};
    flag_select = single_search(target_cost, index, temp_index, search_domain[0].begin(), 
                  search_vector, search_domain, flag, temp_count, total_conut, tempvector,
                  cleanup_queue,
                  min_cost, min_cost_index);
    if(flag_select) break;
  }
  return flag_select;
}

uint MatchingGenerator::print_statics(bool print_flag) {
  uint count[rs_n + 1] = {0};
  bool flag[stripes_num] = {0};
  for (uint i = 0; i < stripes_num; ++i) {
    if (!flag[i]){
      bool a = true;
      for(auto j : stripes[i].matching_index){
        //cout << j << endl;
        if( j > stripes_num) break; 
        if(flag[j]){
          a = false;
          break;
        }
      }
      if(a){
        ++count[stripes[i].cur_cost];
        flag[i] = true;
        for(auto j : stripes[i].matching_index){
          if( j > stripes_num) break;
          flag[j] = true;
        }
      }
    }
  }
  
  if (print_flag) {
    std::cout << "--- zero rate = " << 200.0 * count[0] / stripes_num << "%\n";
  }
  uint total_costs = 0;
  for (uint8_t i = 0; i <= rs_n; ++i) {
    if (!count[i]) {
      continue;
    }
    total_costs += i * count[i];
    if (print_flag) {
      std::cout << "cost = " << +i << " :  " << count[i] << '\n';
    }
  }
  if (print_flag) {
    std::cout << "Total_Costs = " << total_costs << '\n';
  }
  return total_costs;
}

void MatchingGenerator::dist_based_search() {
  if (node_num < 2 * rs_k + rs_m) {
    std::cerr << "bad input parameter!\n";
    exit(-1);
  }
  build_map();
  struct timeval start_time_0, end_time_0, start_time_1, end_time_1;
  //gettimeofday(&start_time_0, nullptr);
  for(uint cluster = 0 ; cluster < K; cluster++){
    cout << cluster <<" : " << parity_cluster_0[cluster].size() << endl;
    for (uint index = 0; index < parity_cluster_0[cluster].size(); ++index){
      remaining_stripes.push_back(parity_cluster_0[cluster][index]);
    }
    std::queue<std::list<uint>::iterator> it_to_del;
    for (uint8_t cost = 0; cost < x * (rs_m - 1) && !remaining_stripes.empty(); ++cost) {
      //cout << "target cost: " << +cost << endl;
      for (auto cur_it = remaining_stripes.begin();
            cur_it != remaining_stripes.end(); ++cur_it) {
        uint &cur_stripe = *cur_it;
        //cout << "cur_stripe: " << cur_stripe << endl;
        if (stripes[cur_stripe].selected) {
          /*
          cout << "cluster: "<< cluster << " " << cur_stripe << " : " ;
          for(auto it :stripes[cur_stripe].matching_index){
            cout << it  << " ";
          }
          cout << endl;
          */
          //it_to_del.push(cur_it);
          continue;
        }
        bool succeed = false;
        //cout << "stripes[cur_stripe].cur_cost: " << +stripes[cur_stripe].cur_cost << " cost " << +cost << endl;
        if (stripes[cur_stripe].cur_cost == cost) {
          succeed = true;
          for(auto i : stripes[cur_stripe].matching_index){
            if(stripes[i].selected){
              succeed = false;
              break;
            }
          }
          if (succeed);
          else {
            stripes[cur_stripe].matching_index = {UINT32_MAX};
            stripes[cur_stripe].cur_cost = UINT8_MAX;
            succeed = search_for_matching(cur_stripe, cost ,cluster);
          }
        } else {
          succeed = search_for_matching(cur_stripe, cost, cluster);
        }
        if (succeed) {
          mark_matched(cur_stripe, stripes[cur_stripe].matching_index,
                        stripes[cur_stripe].cur_cost);
          //it_to_del.push(cur_it);
        }
      }
      for (auto cur_it = remaining_stripes.begin();
            cur_it != remaining_stripes.end(); ++cur_it) {
        uint &cur_stripe = *cur_it;
        if (stripes[cur_stripe].selected)
          it_to_del.push(cur_it);
      }
      while (!it_to_del.empty()) {
        auto x = it_to_del.front();
        it_to_del.pop();
        //cout << *x << " ";
        remaining_stripes.erase(x);
      }
    }
    // gettimeofday(&end_time_0, nullptr);
    // double fin_time_0 = (end_time_0.tv_sec - start_time_0.tv_sec) * 1000 +
    //             (end_time_0.tv_usec - start_time_0.tv_usec) / 1000;
    // cout << "finish time_0 : " << fin_time_0 << "ms\n";
    /*
    for(uint i = 0 ; i < stripes_num ; i++)
    cout << stripes[i].selected << " " ;
    cout << endl;
    for(auto i : remaining_stripes)
      cout << i << " " ;
    cout << endl;
    */
    cout << "remaining_stripe's num = "<< remaining_stripes.size() << endl;
    //gettimeofday(&start_time_1, nullptr);
    if (!remaining_stripes.empty()) {
      basic_greedy(true);
    }
    // gettimeofday(&end_time_1, nullptr);
    // double fin_time_1 = (end_time_1.tv_sec - start_time_1.tv_sec) * 1000 +
    //             (end_time_1.tv_usec - start_time_1.tv_usec) / 1000;
    // cout << "finish time_1 : " << fin_time_1 << "ms\n";
    //cout << remaining_stripes.size() << endl;
  }
}

uint MatchingGenerator::get_distance(std::vector<uint> &points, std::vector<uint> &center){
  std::unordered_set<uint> result_set; // 存放结果
  unordered_set<uint> nums_set(points.begin(), points.end()); 
  for (uint it : center) {
      if (nums_set.find(it) != nums_set.end()) {
          result_set.insert(it);
      }
  }
  return result_set.size();
}

void MatchingGenerator::k_means(bool flag_dist_search) {
  K = 1;
  if(!flag_dist_search){
    for(uint i = 0 ; i < stripes_num ; i++){
      parity_cluster_0.resize(1);
      parity_cluster_0[0].push_back(i);
    }
    return;
  }
  uint p = x;
  while(p--){
    K *= x;
  }
  uint inters = 20;
  center_cluster.resize(K);
  srand(time(NULL));
  for(uint i = 0 ; i < K; i++){
    uint j = rand() % stripes_num;
    center_cluster[i].assign(stripes[j].parity_dist.begin(),stripes[j].parity_dist.end());
  }
  std::vector<uint> label(stripes_num,-1);
  for(uint i = 0 ; i < inters ; i++){
    for(uint pointnum = 0 ; pointnum < stripes_num; pointnum++){
      uint distance = 0;
      srand(time(NULL));
      for(uint cluster = 0; cluster < K; cluster++){
        std::vector<uint> a(rs_m);
        a.assign(stripes[pointnum].parity_dist.begin(),stripes[pointnum].parity_dist.end());
        uint temp_distance = get_distance(a, center_cluster[cluster]);
        if (temp_distance > distance){
					distance = temp_distance;
					label[pointnum] = cluster;
        }
      }
      if(distance == 0){
        label[pointnum] = rand() % K;
      }
    }
    for(uint cluster = 0; cluster < K; cluster++){
			uint count[node_num] = {0};
      vector<vector<uint> > max_node(2, vector<uint>(rs_m, 0));
			for (uint pointnum = 0; pointnum < stripes_num; pointnum++){
				if (label[pointnum] == cluster){ 
          for(uint index = 0; index < rs_m ; index++){
					  ++count[stripes[pointnum].parity_dist[index]];
          }
				}
			}
      for(uint index = 0; index < node_num ; index++){
        uint temp_count = 0;
        while(temp_count < rs_m){
          if(count[index] > max_node[1][temp_count])
            temp_count++;
          else break;
        }
        if(temp_count == 0)
          continue;
        for(uint i = 1; i < temp_count ; i++){
          max_node[0][i-1] = max_node[0][i];
          max_node[1][i-1] = max_node[1][i];
        }
        max_node[0][temp_count-1] = index;
        max_node[1][temp_count-1] = count[index];
      }
			for(uint index = 0; index < rs_m; index++){
        center_cluster[cluster][index] = max_node[0][index];
      }
		}
  }
  /*
  for(auto &a : center_cluster){
    for(auto b : a){
      cout << b << " ";
    }
    cout << endl;
  }
  for(auto i : label){
    cout << i <<" ";
  }
  */

  parity_cluster_0.resize(K);
  for(uint i = 0 ; i < stripes_num ; i++){
    parity_cluster_0[label[i]].push_back(i);
  }
  /*
  for(uint i = 0 ; i < parity_cluster_0.size(); i++){
    cout << " cluster: " << i << " sritpe: ";
    for(auto it : parity_cluster_0[i]){
      cout << " " << it << " ";
    }
    cout << endl;
  }
  */
  if(flag_dist_search == true) return;
}

void MatchingGenerator::parity_means() {
  parity_cluster.resize(node_num);
  for(uint i = 0 ; i < stripes_num; i++){
    for(int index = 0; index < rs_m ; index++){
      if(stripes[i].parity_dist[index]){
        parity_cluster[stripes[i].parity_dist[index]].push_back(i);
      }
    }
  }
}


uint8_t MatchingGenerator::get_cost(std::vector<uint> &merge_stripes_index) {
  /*
  Stripe a = stripes[merge_stripes_index[0]];
  Stripe b = stripes[merge_stripes_index[1]];
  uint8_t cost = 0;
  bool flag[node_num] = {false};
  uint8_t i;
  // check the cost caused by the excess of data chunks
  for (i = 0; i < rs_k; ++i) {
    if (flag[a.data_dist[i]]) {
      ++cost;
    } else {
      flag[a.data_dist[i]] = true;
    }
    if (flag[b.data_dist[i]]) {
      ++cost;
    } else {
      flag[b.data_dist[i]] = true;
    }
  }
  // check the cost of misaligned parity blocks
  //     & the overflow after determining the new parity blocks' placement
  char16_t parity_dest;
  for (i = 0; i < rs_m; ++i) {
    if (a.parity_dist[i] == b.parity_dist[i]) {
      parity_dest = a.parity_dist[i];
    } else {
      ++cost;
      parity_dest =
          flag[a.parity_dist[i]] ? b.parity_dist[i] : a.parity_dist[i];
    }
    if (flag[parity_dest]) {
      ++cost;
    } else {
      flag[parity_dest] = true;
    }
  }
  return cost;
  */
  bool flag[node_num] = {false};
  //std::vector<uint8_t> cost(3,0);
  uint8_t i;
  uint8_t temp_cost = 0;
  
  // check the cost caused by the excess of data chunks
  for (int i = 0; i < rs_k; ++i) {
    for(int k = 0; k < merge_stripes_index.size(); k++){
      if (flag[stripes[merge_stripes_index[k]].data_dist[i]]) {
        ++temp_cost;
      } else {
        flag[stripes[merge_stripes_index[k]].data_dist[i]] = true;
      }
    }
  }
  //cost[1] = temp_cost;
  
  //std::vector<std::vector<bool>> parity_location(rs_m, vector<bool> (node_num, false));

  // check the cost of misaligned parity blocks
  //     & the overflow after determining the new parity blocks' placement
  
  for (i = 0; i < rs_m; ++i) {
    bool parity_flag = true;
    bool parity_location[node_num] = {false};
    for(uint k =0; k < merge_stripes_index.size(); ++k){
      const std::u16string &pd = stripes[merge_stripes_index[k]].parity_dist;
      if(parity_location[pd[i]] == false){
        parity_location[pd[i]] = true;
        temp_cost++;
        if(flag[pd[i]] == false)
          parity_flag = false;
      }
    }
    if(!parity_flag)  temp_cost--;
  }
  
  //cost[0] = temp_cost;
  //cost[2] = temp_cost - cost[1];
  return temp_cost;
  /*
  for (uint k = 0; k < merge_stripes_index.size(); ++k) {
    const std::u16string &pd = stripes[merge_stripes_index[k]].parity_dist;
    for (uint i = 0 ; i < pd.length(); i++){
      parity_location[i][pd[i]] = true;
    }
  }
  for(int k = 0; k < rs_m; k++){
    bool parity_flag = true;
    uint8_t count = 0;
    for(int n = 0; n < node_num; n++){
      if(parity_location[k][n] == true){
          count++;
          if(flag[n] == false)
            parity_flag = false;
      }
    }
    int count1 = count;
    parity_flag ? temp_cost += count : temp_cost += count - 1;
  }
  cost[0] = temp_cost;
  cost[2] = temp_cost - cost[1];
  return cost;
  */
}

void MatchingGenerator::input_stripes(char16_t **raw, char16_t **parity) {
  //cout << "stripes_num:" << " " << stripes_num << endl;
  for (uint i = 0; i < stripes_num; ++i) {
    //cout << "stripe " << i <<": ";
    stripes[i].set_info(i, raw[i], parity[i]);
    /*
    for(int j = 0; j < rs_k; j++)
      cout << raw[i][j] << " ";
    for(int j = 0; j < rs_m; j++)
      cout << parity[i][j] << " ";
    cout<< endl;
    */
  }
}

uint MatchingGenerator::basic_greedy(bool remain_flag) {
  if (node_num < x * rs_k + rs_m) {
    std::cerr << "bad input parameter!\n";
    exit(-1);
  }
  //struct timeval start_time, end_time, start_time_0, end_time_0;
  uint8_t cost;
  uint i, j, a, b;
  std::vector<std::vector<uint>> level[rs_n + 1];
  std::vector<std::vector<std::list<uint>::iterator>> level_dist[rs_n + 1];
  double fin_time_0 =  0.0;
  double fin_time_2 =  0.0;
  double fin_time_3 =  0.0;
  if(!remain_flag){
    std::vector<uint> it_vector;
    uint i = 0;
    it_vector.emplace_back(i);
    int total = 0;
    while(i < stripes_num){
      if(it_vector.size() == x){
        cost = get_cost(it_vector);
        level[cost].emplace_back(it_vector);
        total++;
        i = it_vector.back();
        it_vector.pop_back();
      }
      else{
        i++;
        if(i < stripes_num){
          it_vector.push_back(i);
        }
        else{
          if(it_vector.size()){
            i = it_vector.back();
            it_vector.pop_back();
          }
        }
      }
    }
    
    // gettimeofday(&end_time_0, nullptr);
    // double fin_time_1 = (end_time_0.tv_sec - start_time_0.tv_sec) * 1000 +
    //             (end_time_0.tv_usec - start_time_0.tv_usec) / 1000;
    // cout << "finish time _1: " << fin_time_1 << "ms\n";
    //cout << "finish time _2: " << fin_time_2 << "ms\n";
    //cout << "finish time _3: " << fin_time_3 << "ms\n";
  }
  else{
    std::vector<std::list<uint>::iterator> it_vector;
    std::list<uint>::iterator it = remaining_stripes.begin();
    bool flag_it[stripes_num] = {0};
    it_vector.emplace_back(it);
    flag_it[*it] = true;
    uint k = 1;
    while(it != remaining_stripes.end()){
      if(k == x){
        std::vector<uint> tempstripe;
        for(auto it_0 : it_vector){
          tempstripe.emplace_back(*it_0);
        }
        cost = get_cost(tempstripe);
        level_dist[cost].emplace_back(it_vector);
        it = it_vector.back();
        it_vector.pop_back();
        flag_it[*it] = false;
        k--;
      }
      else{
        it++;
        if(it != remaining_stripes.end()){
          k++;
          it_vector.emplace_back(it);
          flag_it[*it] = true;
        }
        else{
          if(it_vector.size()){
            it = it_vector.back();
            it_vector.pop_back();
            flag_it[*it] = false;
            k--;
          }
        }
      }
    }
  }
  uint selected[rs_n + 1] = {0};
  uint num_to_select;
  num_to_select =
      remain_flag ? remaining_stripes.size() / x : stripes_num / x;
  bool flag[stripes_num] = {0};
  if(!remain_flag){
    for (uint8_t cur_cost = 0; cur_cost <= rs_n && num_to_select; ++cur_cost) {
      for (auto &m : level[cur_cost]){
        bool total_flag = true;
        for(auto it : m){
          if(!flag[it]) continue;
          else {
            total_flag = false;
            break;
          }
        }
        if(total_flag){
          ++selected[cur_cost];
          for(auto it : m){
            flag[it] = true;
          }
          std::vector<uint> temp_vector;
          //i = m[0];
          //temp_vector.assign(m.begin() + 1, m.end());
          //mark_matched(i, temp_vector, cur_cost, false);
          if (!(--num_to_select)) {
            break;
          }
        }
      }
    }
  }
  else{
    for (uint8_t cur_cost = 0; cur_cost <= rs_n && num_to_select; ++cur_cost) {
      if(!level_dist[cur_cost].size()) continue;
      for (auto &m : level_dist[cur_cost]) {
        i = *m[0];
        j = *m[1];
        if(2 < x){
          a = *m[2];
          if(3 < x){
            b = *m[3];
          }
        }
        
        if(x == 2){
          if (!flag[i] && !flag[j]){
            ++selected[cur_cost];
            flag[i] = flag[j] = true;
            std::vector<uint> temp_vector = {j};
            mark_matched(i, temp_vector, cur_cost, false);
            remaining_stripes.erase(m[0]);
            remaining_stripes.erase(m[1]);
            if (!(--num_to_select)) {
              break;
            }
          }
        }
        if(x == 3){
          if (!flag[i] && !flag[j] && !flag[a]){
            ++selected[cur_cost];
            flag[i] = flag[j] = flag[a] = true;
            std::vector<uint> temp_vector = {j,a};
            mark_matched(i, temp_vector, cur_cost, false);
            remaining_stripes.erase(m[0]);
            remaining_stripes.erase(m[1]);
            remaining_stripes.erase(m[2]);
            if (!(--num_to_select)) {
              break;
            }
          }
        }
        if(x == 4){
          if (!flag[i] && !flag[j] && !flag[a] && !flag[b]){
            ++selected[cur_cost];
            std::vector<uint> temp_vector = {j,a,b};
            flag[i] = flag[j] = flag[a] = flag[b]= true;
            mark_matched(i, temp_vector, cur_cost, false);
            remaining_stripes.erase(m[0]);
            remaining_stripes.erase(m[1]);
            remaining_stripes.erase(m[2]);
            remaining_stripes.erase(m[3]);
            if (!(--num_to_select)) {
              break;
            }
          }
        }
      }
    }
  }
  if (!remain_flag) {
    uint total_costs = 0;
    for (uint8_t cur_cost = 0; cur_cost <= rs_n; ++cur_cost) {
      if (!selected[cur_cost]) {
        continue;
      }
      total_costs += cur_cost * selected[cur_cost];
      cout << "cost = " << +cur_cost << " : " << selected[cur_cost] << endl;
    }
    return total_costs;
  } 
  else {
    return 0;
  }
}
void MatchingGenerator::random_match() {
  std::vector<uint> pos;
  uint cost_count_table[rs_n + 1] = {0};
  for (uint i = 0; i < stripes_num; ++i) {
    pos.push_back(i);
  }
  //srand(time(nullptr));
  uint cost = 0;
  std::random_shuffle(pos.begin(), pos.end());
  for (uint i = 0; i < stripes_num; i += x) {
    uint8_t temp_cost = 0;
    std::vector<uint> tempstripe;
    for(uint j = 0; j < x ; j++){
      tempstripe.push_back(pos[i+j]);
    }
    /*
    for(auto it : tempstripe){
      cout << it <<" ";
    }
    cout << endl;
    */  
    ++cost_count_table[get_cost(tempstripe)];
    tempstripe.clear();
  }
  for(uint8_t i = 1; i < rs_n + 1; i++){
    cost += i * cost_count_table[i];
  }
  cout << "cost:" << cost << endl;
}