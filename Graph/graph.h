#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <stack>
#include <numeric> // for accumulate
#include <string>
#include <dirent.h>
#include <climits>
#include <cmath>
#include <omp.h>
#include <map>
#include <iostream>

using namespace std;

// for reading graph files
int read_graphs(string &fname, int &num_nodes, int &num_edges, vector<int> &row_ptr, vector<int> &col_ind);

// for getting file names in a directory
void get_filenames(vector<string> &filenames, const vector<string> &locations);

// to perform all the computations at once
void compute_centralities(const int &num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, 
                          vector<vector<float>> &result, const vector<bool> &to_calculate, int step_size=1000);
void comp_cent(const int &num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, 
                          vector<vector<float>> &result, const vector<bool> &to_calculate, int step_size=1000);
void print_centralities(string filename, vector<vector<float>> &result);

// bfs algorithm sequential implementation
void bfs(int start_node, int num_nodes, vector<int> row_ptr, vector<int> col_ind, vector<int> &distance_arr, int step_size = INT_MAX);

// centrality functions
void closeness_centrality_approx(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering, int size = 1000);
void closeness_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering);
void degree_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering);
void degree2_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering);

#endif
