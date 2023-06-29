#ifndef STRIPE_HH
#define STRIPE_HH

#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <unordered_set>

class Stripe {
  static uint16_t node_num;
  static uint8_t rs_k;
  static uint8_t rs_m;

 public:
  uint stripe_index;
  std::vector<uint> matching_index;
  uint8_t cur_cost;
  std::u16string parity_dist;
  std::u16string data_dist;
  bool selected;

  Stripe();
  ~Stripe();
  void set_info(uint _stripe_idnex, char16_t *data_dist_src,
                char16_t *parity_dist_src);
  static void set_params(uint16_t _node_num, uint8_t k, uint8_t m);
  static void generate_random_distributions(char16_t **raw_dist,
                                            char16_t **parity_dist,
                                            uint16_t node_num, uint stripes_num,
                                            uint8_t rs_k, uint8_t rs_m);
};

#endif