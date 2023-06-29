#include <sys/time.h>

#include "matching_generator.hh"

using namespace std;

#define COMPARE true

int main(int argc, char **argv) {
  if (argc != 6) {
    cerr << "Usage: ./matching_main x [stripes_num] [node_num] [rs_k] [rs_m]"
         << endl;
    exit(-1);
  }
  uint x = strtoul(argv[1], nullptr, 10);
  uint stripes_num = strtoul(argv[2], nullptr, 10);
  uint16_t node_num = strtoul(argv[3], nullptr, 10);
  uint8_t rs_k = strtoul(argv[4], nullptr, 10);
  uint8_t rs_m = strtoul(argv[5], nullptr, 10);
  char16_t **raw_dist = new char16_t *[stripes_num];
  char16_t **parity_dist = new char16_t *[stripes_num];
  for (uint i = 0; i < stripes_num; ++i) {
    raw_dist[i] = new char16_t[rs_k];
    parity_dist[i] = new char16_t[rs_m];
  }
  Stripe::generate_random_distributions(raw_dist, parity_dist, node_num,
                                        stripes_num, rs_k, rs_m);
  
  /*
  for(uint index = 0; index < stripes_num; index++){
    cout << "stripe " << index << " ";
    for(uint8_t i = 0; i < rs_k ; i++){
      cout << raw_dist[index][i] << " ";
    }
    for(uint8_t j = 0; j < rs_m ; j++){
      cout << parity_dist[index][j] << " ";
    }
    cout << endl;
  }
  */
  struct timeval start_time, end_time;
  double fin_time[6] = {0.0};
  
  cout << "=========== random_match ===========\n";
  gettimeofday(&start_time, nullptr);
  MatchingGenerator *mg2 =
      new MatchingGenerator(x, node_num, rs_k, rs_m, stripes_num);
  mg2->input_stripes(raw_dist, parity_dist);
  mg2->random_match();
  gettimeofday(&end_time, nullptr);
  delete mg2;
  fin_time[1] = (end_time.tv_sec - start_time.tv_sec) * 1000 +
                (end_time.tv_usec - start_time.tv_usec) / 1000;
  cout << "finish time: " << fin_time[1] << "ms\n";
  
  if (COMPARE) {
  cout << "=========== xStripeMerge-G ===========\n";
  gettimeofday(&start_time, nullptr);
  MatchingGenerator *mg3 =
      new MatchingGenerator(x, node_num, rs_k, rs_m, stripes_num);
  mg3->input_stripes(raw_dist, parity_dist);
  uint tot_costs = mg3->basic_greedy();
  gettimeofday(&end_time, nullptr);
  delete mg3;
  fin_time[2] = (end_time.tv_sec - start_time.tv_sec) * 1000 +
                (end_time.tv_usec - start_time.tv_usec) / 1000;
  cout << "finish time: " << fin_time[2] << "ms\n";
  cout << "Total_cost:" << tot_costs << endl;
  cout << "Average_cost:" << double(tot_costs)/double(stripes_num)*x << endl;
  // cout << "\n*** Speedup ratio = " << time_2 / time_1 << " ***\n";
  }
  
  cout << "===========  xStripeMerge-P ===========\n";
  gettimeofday(&start_time, nullptr);
  MatchingGenerator *mg6 =
      new MatchingGenerator(x, node_num, rs_k, rs_m, stripes_num);
  mg6->input_stripes(raw_dist, parity_dist);
  mg6->k_means(false);
  mg6->dist_based_search();
  gettimeofday(&end_time, nullptr);
  uint cost_5 = mg6->print_statics();
  cout << "Average_cost:" << double(cost_5)/double(stripes_num)*x << endl;
  fin_time[5] = (end_time.tv_sec - start_time.tv_sec) * 1000 +
                (end_time.tv_usec - start_time.tv_usec) / 1000;
  cout << "finish time: " << fin_time[5] << "ms\n";
  delete mg6;
  /*
  cout << "=========== xRSmerge_k_means_dist_search  ===========\n";
  gettimeofday(&start_time, nullptr);
  MatchingGenerator *mg7 =
      new MatchingGenerator(x, node_num, rs_k, rs_m, stripes_num);
  mg7->input_stripes(raw_dist, parity_dist);
  mg7->k_means(true);
  mg7->dist_based_search();
  gettimeofday(&end_time, nullptr);
  uint cost_6 = mg7->print_statics();
  fin_time[3] = (end_time.tv_sec - start_time.tv_sec) * 1000 +
                (end_time.tv_usec - start_time.tv_usec) / 1000;
  cout << "finish time: " << fin_time[3] << "ms\n";
  delete mg7;
  */
  for (uint i = 0; i < stripes_num; ++i) {
    delete[] raw_dist[i];
    delete[] parity_dist[i];
  }
  delete[] raw_dist;
  delete[] parity_dist;

  return 0;
}
