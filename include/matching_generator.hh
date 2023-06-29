#ifndef MATCHING_GENERATOR_HH
#define MATCHING_GENERATOR_HH

#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <unordered_map>
#include <set>
#include <vector>
#include <stack>
#include "stripe.hh"

#define UINT32_MASK 0xFFFFFFFF

class MatchingGenerator {
  uint x;
  uint16_t node_num;
  uint8_t rs_k;
  uint8_t rs_m;
  uint8_t rs_n;
  uint stripes_num;

  std::unordered_map<uint, std::set<std::vector<int>>> table;

  uint K ; //x * node_num / rs_m;
  std::vector<Stripe> stripes;
  std::vector<uint8_t> cost_table;
  std::vector<std::vector<uint>> center_cluster;
  std::list<uint> remaining_stripes;
  std::list<uint> final_remaining_stripes;
  std::vector<std::unordered_map<std::u16string, std::list<uint>>> parity_map;
  std::vector<std::unordered_map<std::u16string, std::list<uint>>> partial_dist_map;

  std::vector<std::vector<uint>> parity_cluster;
  std::vector<std::vector<uint>> parity_cluster_0;

  std::unordered_map<std::u16string, std::list<uint>> zero_map;
  std::vector<uint> none_zero_candidate;

 public:
  MatchingGenerator() {}
  ~MatchingGenerator() {}
  MatchingGenerator( uint8_t stripe_x, uint16_t _node_num, uint8_t k, uint8_t m,
                    uint _stripes_num);
  void input_stripes(char16_t **raw, char16_t **parity);
  uint8_t get_cost(std::vector<uint> &merge_stripes_index);
  uint8_t get_total_cost();
  uint get_distance(std::vector<uint> &points, std::vector<uint> &center);
  uint print_statics(bool print_flag = true);
  uint basic_greedy(bool remain_flag = false);
  void graph_based_search();

  void random_match();
//  long long blossom();
  void parity_means();
  void k_means(bool flag_dist_search);
  void build_map_test();
  void dist_based_search();

 private:
  std::unordered_map<uint, std::set<std::vector<int>>> build_search_table(int i);
  void build_map();
  std::vector<uint8_t> get_cost(uint64_t token);
  void mark_matched(uint i, std::vector<uint> j, uint8_t cost, bool flag = true);
  bool single_search(uint8_t target_cost, uint index, uint temp_index, std::list<uint>::iterator temp_it,  std::vector<int> search_vector, 
                    std::vector<std::list <uint>> search_domain, std::vector<bool> flag, std::vector<int> temp_count,
                    uint total_conut,
                    std::vector<uint> tempvector,
                    std::queue<uint> cleanup_queue,
                    uint8_t min_cost,
                    std::vector<uint> min_cost_index);
  bool search_for_matching(uint index, uint8_t target_cost, uint cluseter);
//  void make_optimization();
//  void get_one_scheme(uint index);
};

#endif
