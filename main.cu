/*
 @author furkanergun
 @author barisbatuhan
 @author halpar
 @author caghankoksal
                        @sabanciuniv.edu
 */
#include <iostream>
#include <cuda.h>
#ifndef UGRAPH_H
#define UGRAPH_H

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <dirent.h>
#include <cmath>
#include <omp.h>
#include <iostream>

class UGraph
{
public:
    
    std::vector<int> row_ptr;
    std::vector<int> col_ind;
    int num_nodes;
    int num_edges;
    int * rowPtr;
    int * colPtr;
    // getter
    int getNumNodes(){
        return num_nodes;
    }
    
    int getNumEdges(){
        return num_edges;
    }
    int * getRowPtr(){
        //allocates space
        int * rowPtr = new int[row_ptr.size()];
        int * itr = rowPtr;
        for(auto &i : row_ptr){
            *itr=i;
            itr++;
        }
        return rowPtr;
    }
    int * getColPtr(){
        //allocates space
        int * colPtr = new int[col_ind.size()];
        int * itr = colPtr;
        for(auto &i : col_ind){
            *itr=i;
            colPtr++;
        }
        return colPtr;
    }
    // constructor
    UGraph(std::string fname);
    UGraph(int node_cnt, int edge_cnt); // random graph generation

    // centrality function
    void compute_centralities(std::vector<std::vector<float>> &result, const std::vector<bool> &requested);
    
    // centrality helper functions
    int bfs_topdown(std::vector<int> &dist, int source, std::vector<double> &sigma, std::vector<int> &queue,
                    std::vector<int> &dist_counter, const std::vector<bool> &requested, int &queue_size);
    float node_closeness(std::vector<int> &dist_counter, int max_dist);
    void node_betweenness(std::vector<std::vector<float>> &result, std::vector<int> &queue,
                          std::vector<int> &dist, std::vector<double> &sigma, int &queue_size);

    // static methods
    static void get_filenames(std::vector<std::string> &filenames, const std::vector<std::string> &locations);
    static void print_centralities(std::string filename, std::vector<std::vector<float>> &result);

private:
    std::string family;
    std::string relative_path;
};

/* CONSTRUCTOR METHODS */

UGraph::UGraph(std::string fname)
{
    std::ifstream input(fname.c_str());
    if (input.fail())
    {
        throw "No file is found in the current path!";
    }
    else
    {
        relative_path = fname;
    }
    // read Ugraph
    std::string line = "%";
    family = "%";
    while (line.find("%") != std::string::npos)
    {
        getline(input, line);
        if (family == "%" && line.find("kind:") != std::string::npos)
        {
            family = line.substr(8);
            family = family.substr(0, family.length() - 1);
        }
    }

    std::istringstream ss(line);
    ss >> num_nodes >> num_nodes >> num_edges;
    int v1, v2;
    double weight;

    std::vector<std::vector<int>> adj_list(num_nodes);
    for (int i = 0; i < num_edges; i++)
    {
        getline(input, line);
        std::istringstream inp(line);
        inp >> v1 >> v2;
        v1--; // make it 0 based
        v2--;

        if (v1 != v2)
        {
            adj_list[v1].push_back(v2); // add the edge v1->v2
            adj_list[v2].push_back(v1); // add the edge v2->v1
        }
    }

    row_ptr = std::vector<int>(num_nodes + 1);
    col_ind = std::vector<int>(2 * num_edges);
    row_ptr[0] = 0;
    int index = 0;
    for (int v = 0; v < num_nodes; v++)
    {
        row_ptr[v + 1] = adj_list[v].size(); // assign number of edges going from node v
        for (int i = 0; i < (int)adj_list[v].size(); i++)
        {
            col_ind[index] = adj_list[v][i]; // put all edges in order wrt row_ptr
            index++;
        }
    }
    for (int v = 1; v < num_nodes + 1; v++)
    { // cumulative sum
        row_ptr[v] += row_ptr[v - 1];
    }
}

UGraph::UGraph(int node_cnt, int edge_cnt)
{
    srand(112);
    std::vector<std::vector<bool>> adj_list(node_cnt, std::vector<bool>(node_cnt, false));
    num_nodes = node_cnt;
    num_edges = edge_cnt;
    int edge_num = edge_cnt;
    while (edge_num > 0)
    {
        int node1 = rand() % node_cnt;
        int node2 = rand() % node_cnt;
        if (node1 == node2)
            continue;
        if (adj_list[node1][node2] == false)
        {
            adj_list[node1][node2] = true;
            adj_list[node2][node1] = true;
            edge_num--;
        }
    }
    row_ptr = std::vector<int>(node_cnt + 1);
    col_ind = std::vector<int>(2 * edge_cnt);
    row_ptr[0] = 0;
    int index = 0;
    for (int v = 0; v < node_cnt; v++)
    {
        int adj_cnt = 0;
        for (int i = 0; i < (int)adj_list[v].size(); i++)
        {
            if (adj_list[v][i] == true)
            {
                col_ind[index] = i; // put all edges in order wrt row_ptr
                index++;
                adj_cnt++;
            }
        }
        row_ptr[v + 1] = row_ptr[v] + adj_cnt; // assign number of edges going from node v
    }
}

/* CENTRALITY METHODS */

int UGraph::bfs_topdown(std::vector<int> &dist, int source, std::vector<double> &sigma, std::vector<int> &queue,
                        std::vector<int> &dist_counter, const std::vector<bool> &requested, int &queue_size)
{
    queue[0] = source;
    dist[source] = 0;
    sigma[source] = (double) 1 / (num_nodes * num_nodes);
    dist_counter[0]++;
    int front = 0, max_distance = 0;
    queue_size = 1;

    while(front < queue_size) {
        int &v = queue[front];
        if(dist[v] == 2 && !requested[2] && !requested[3]) break;
        front++;
        for(int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
            int &w = col_ind[edge];
            if(dist[w] < 0) { // w not visited yet
                queue[queue_size] = w;
                queue_size++;
                dist[w] = dist[v] + 1;
                if(max_distance < dist[w]) max_distance = dist[w];
                dist_counter[dist[w]]++;
            }
            if(requested[3] && dist[w] == dist[v] + 1) {
                sigma[w] += sigma[v];
            }
        }
    }

    return max_distance;
}

float UGraph::node_closeness(std::vector<int> &dist_counter, int max_dist)
{
    int sum = 0;
    for(int i = 1; i <= max_dist; i++) {
        sum += i * dist_counter[i];
    }
    return (float) 1 / sum;
}

void UGraph::node_betweenness(std::vector<std::vector<float>> &result, std::vector<int> &queue,
                              std::vector<int> &dist, std::vector<double> &sigma, int &queue_size)
{
    std::vector<float> delta(num_nodes, 0);
    for(int index = queue_size - 1; index > 0; index--) {
        int &w = queue[index];
        if(w == -1) continue;
        for(int edge = row_ptr[w]; edge < row_ptr[w + 1]; edge++) {
            int &v = col_ind[edge];
            if(dist[v] + 1 == dist[w]) { // v is parent of w
                delta[v] += (float) (sigma[v] / sigma[w]) * (1 + delta[w]);
            }
        }
        #pragma omp atomic
        result[3][w] += delta[w];
    }
    return;
}

void UGraph::compute_centralities(std::vector<std::vector<float>> &result, const std::vector<bool> &requested)
{
    result = std::vector<std::vector<float>>(4, std::vector<float>(num_nodes, 0));

    #pragma omp parallel for shared(result) schedule(dynamic)
    for(int s = 0; s < num_nodes; s++) { // for all source vertices
        // Parameter Initialization
        std::vector<int> dist(num_nodes, -1);            // distance array from source to vertices
        std::vector<int> dist_counter(num_nodes, 0);     // counts nodes for each distances
        std::vector<double> sigma(num_nodes, 0);            // number of shortest paths from source to index
        std::vector<int> queue(num_nodes, -1);           // non-decreasing read order for bfs
        int queue_size = 0;
        
        // Degree 1 calculation
        if(requested[0]) result[0][s] = (float) row_ptr[s + 1] - row_ptr[s];
        // BFS - maximum distance (diameter of graph) is returned
        int max_dist = bfs_topdown(dist, s, sigma, queue, dist_counter, requested, queue_size);
        // Degree 2 calculation
        if(requested[1]) result[1][s] = (float) dist_counter[2];
        // Closeness calculation
        if(requested[2]) result[2][s] = node_closeness(dist_counter, max_dist);
        // Betweenness calculation
        if(requested[3]) node_betweenness(result, queue, dist, sigma, queue_size);
    }
}

/* STATIC METHODS */
void UGraph::get_filenames(std::vector<std::string> &filenames, const std::vector<std::string> &locations)
{
    for (int i = 0; i < locations.size(); i++)
    {
        if (auto dir = opendir(locations[i].c_str()))
        {
            while (auto f = readdir(dir))
            {
                if (!f->d_name || f->d_name[0] == '.')
                    continue;

                std::string path = locations[i] + f->d_name;
                filenames.push_back(path);
            }
        }
    }
}

void UGraph::print_centralities(std::string filename, std::vector<std::vector<float>> &result) {
    std::cout << std::endl
         << "--------------------------------------------------------------------" << std::endl
         << "| Results for " << filename << std::endl
         << "--------------------------------------------------------------------" << std::endl
         << "|   Indx     |    Deg1    |    Deg2    |    Clos    |    Betw    |" << std::endl
         << "   ------        ------       ------       ------       ------   " << std::endl;

    for(int i = 0; i < result[0].size(); i++) {
        std::cout << "|";
        std::string index = std::to_string((float) i);
        if(i >= 100) std::cout << " " <<  index << " |";
        else if(i >= 10) std::cout << " " <<  index << "  |";
        else if(i < 0) std::cout << " " <<  index << "  |";
        else std::cout << "  " <<  index << "  |";
        
        for(int j = 0; j < 4; j++) {
            std::string res = std::to_string(result[j][i]);
            // float res = result[j][i];
            if(result[j][i] >= 100) std::cout << " " <<  res << " |";
            else if(result[j][i] >= 10) std::cout << " " <<  res << "  |";
            else if(result[j][i] < 0) std::cout << " " <<  res << "  |";
            else std::cout << "  " <<  res << "  |";
        }
        std::cout << std::endl;
    }
}

#endif

#define mem_per_block 49152
#define mem_const 65536
#define warp_size 32
//max global memory is 11441mb
#define global_mem 11441e6
#define max_thread_per_block 1024
#define INF 99999
#define DEBUG
inline cudaError_t checkCuda(cudaError_t result) {
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
  }
  #endif
  return result;
}

/*
 
 adjecancyListPointers == row_ptr
 adjecancyList == col_ptr
 */

__global__ void betweenKernel(double * result, double * preq, UGraph * graph, int * sigma, int * dist, int *pre, int *prepointer, int no_nodes){
    int tid = threadIdx.x;
    int bDim = blockDim.x;
    if(tid<no_nodes){
        __shared__ int preLen;
        __shared__ int prePtrLen;
        __shared__ int iter;
        if(!tid)
            iter=-1; //init
        __syncthreads();
        while(iter < (no_nodes-1)){
            if(!tid){
                pre[0] = ++iter;
                prePtrLen = 1;
                prePtrLen = 1;
                prepointer[0] = 0;
                prepointer[1] = 1;
            }
            __syncthreads();
            for(int vertex = tid; vertex<no_nodes; vertex+=bDim){
                preq[vertex] = 0.0;
                if(vertex == iter){
                    //same
                    dist[vertex] = 0;
                    sigma[vertex] = (double) 1 / (no_nodes * no_nodes);
                }
                else{
                    dist[vertex] = INF;
                    sigma[vertex] = 0;
                }
            }
            __syncthreads();
            
            bool cont=true;
            while(cont){
                __syncthreads();
                for(int i=tid; i<prepointer[prePtrLen]; i+=bDim){
                    if(i>=prepointer[prePtrLen-1]){
                        int vertex = pre[i]; //no bad access
                        int * rowPtrLoc = graph->rowPtr;
                        for(int j = rowPtrLoc[vertex]; j<rowPtrLoc[vertex+1]; j++){
                            int ed = graph->colPtr[j];
                            int el = atomicCAS(&dist[ed], INF, dist[vertex]+1);
                            if(el == INF){
                                int res = atomicAdd(&preLen, 1);
                                pre[res] = ed;
                            }
                            if((dist[ed]+1) == dist[vertex]){
                                atomicAdd(&sigma[ed], sigma[vertex]);
                            }
                        }
                    }
                }
                __syncthreads();
                if(preLen == prepointer[prePtrLen])
                    cont=false; // no upd
                if(cont && !tid){
                    prePtrLen+=1;
                    prepointer[prePtrLen] = preLen;
                }
                if(cont)
                    __syncthreads();
            }
            
            while(prePtrLen){
                for(int i=tid; i<prepointer[prePtrLen]; i+=bDim){
                    if(i>=prepointer[prePtrLen-1]){
                        int vertex = pre[i];
                        int * rowPtrLoc = graph->rowPtr;
                        for(int j=rowPtrLoc[vertex]; j<rowPtrLoc[vertex+1]; j++){
                            int ed = graph->colPtr[j];
                            if((dist[ed]-1)==dist[vertex])
                                if(sigma[ed])
                                    preq[vertex] += (preq[ed]+1) * (double(sigma[vertex])/sigma[ed]);
                            
                            if(vertex!=iter){
                                result[vertex] += preq[vertex]*(1/2); //adding half
                            }
                                
                        }
                    }
                }
                __syncthreads();
                
                if(!tid)
                    prePtrLen-=1;
                __syncthreads();
            }
        }
    }
}



int main(){
    //host
    UGraph * hostGraph = new UGraph("data/small/epb1.mtx");
    int no_nodes = hostGraph->getNumNodes();
    int no_edges = hostGraph->getNumEdges();
    int * row_ptr = hostGraph->getRowPtr();
    int * col_ptr = hostGraph->getColPtr();
    int rowSize = no_nodes+1;
    int colSize = 2*no_edges+1;
    
    //device
    UGraph * devGraph;
    int devNodes;
    int devEdges;
    int * devRow;
    int * devCol;
    
    
    //centrality host
    double * result = new double[no_nodes];
    
    //centrality device
    int *sigma;
    double *devResult;
    double *preq;
    int *dist;
    int *pre;
    int *prepointer;
    printf("No nodes: %d\n", no_nodes);
    printf("No edges: %d\n", no_edges);
    printf("Started cuda allocations\n");
    cudaMalloc((void **)&devGraph, sizeof(UGraph));
    cudaMalloc((void **)&row_ptr, sizeof(int) * rowSize);
    cudaMalloc((void **)&col_ptr, sizeof(int) * colSize);
    cudaMalloc((void **)&devResult , sizeof(double) * no_nodes);
    cudaMalloc((void **)&preq , sizeof(double) * no_nodes);
    cudaMalloc((void **)&sigma , sizeof(int) * no_nodes);
    cudaMalloc((void **)&dist , sizeof(int) * no_nodes);
    cudaMalloc((void **)&pre , sizeof(int) * (no_nodes+1)); //this two needs additional memory
    cudaMalloc((void **)&prepointer , sizeof(int) * (no_nodes+1));  //this two needs additional memory
    printf("Ended cuda allocations\n");
    printf("Started cuda memory copying\n");
    cudaMemcpy(devGraph, hostGraph, sizeof(UGraph), cudaMemcpyHostToDevice);
    cudaMemcpy(devRow, row_ptr, sizeof(int) * rowSize, cudaMemcpyHostToDevice);
    cudaMemcpy(&(devGraph->rowPtr), &devRow, sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(devCol, col_ptr, sizeof(int) * colSize, cudaMemcpyHostToDevice);
    cudaMemcpy(&(devGraph->colPtr), &devCol, sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(devResult, result, sizeof(double) * no_nodes, cudaMemcpyHostToDevice);
    printf("Ended cuda memory copying\n");
    printf("Starting kernel\n");
    betweenKernel<<<1, 256>>>(devResult, preq, devGraph, sigma, dist, pre, prepointer, no_nodes);
    printf("Kernel ended\n");
    printf("Copying back the result to memory\n");
    cudaMemcpy(result, devResult, sizeof(double) * no_nodes, cudaMemcpyDeviceToHost);
    for(int i=0; i<no_nodes; i++){
        printf("%d: %d\n", i, result[i]);
    }
    printf("Ended");
    cudaFree(devGraph);
    cudaFree(devRow);
    cudaFree(devCol);
    cudaFree(sigma);
    cudaFree(devResult);
    cudaFree(preq);
    cudaFree(dist);
    cudaFree(pre);
    cudaFree(prepointer);
}
