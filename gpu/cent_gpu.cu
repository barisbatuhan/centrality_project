#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cuda.h>

#define THREAD_COUNT 1024
// Max device memory : 4 GB
#define MAX_MEMORY ((long long)4e9)


void read_graph(std::string fname, int *&row_ptr, int *&col_ind, int &num_nodes, int &num_edges, bool zero_based = false)
{
    std::ifstream input(fname.c_str());
    if (input.fail())
        throw "No file is found in the current path!";

    // read graph
    std::string line = "%";
    while (line.find("%") != std::string::npos)
    {
        getline(input, line);
    }

    std::istringstream ss(line.c_str());
    ss >> num_nodes >> num_nodes >> num_edges;
    int edge_cnt = 0;
    int v1, v2;
    std::vector< std::vector<int> > adj_list(num_nodes);
    for (int i = 0; i < num_edges; i++)
    {
        getline(input, line);
        std::istringstream inp(line.c_str());
        inp >> v1 >> v2;
        if (!zero_based)
        {
            v1--; // make it 0 based
            v2--;
        }
        if (v1 != v2)
        {
            adj_list[v1].push_back(v2); // add the edge v1->v2
            adj_list[v2].push_back(v1); // add the edge v2->v1
            edge_cnt++;
        }
    }
    input.close();
    num_edges = edge_cnt;

    cudaMallocHost((void **)&row_ptr, sizeof(int) * (num_nodes + 1));
    cudaMallocHost((void **)&col_ind, sizeof(int) * (2 * num_edges));

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

__global__
void cent_kernel(float *results, int *dist, int *sigma, float *delta, int *rp, int *ci, int n) {
    
    __shared__ int level;
    __shared__ int improved;
    for(int s = blockIdx.x; s < n; s += gridDim.x) {
        if(threadIdx.x == 0) {
	        results[s] = rp[s + 1] - rp[s]; // degree 1
	        level = 0;
	        improved = 1;
                dist[s * n + s] = 0;
                sigma[s * n + s] = 1;
	    }
	    __syncthreads();
        
	    // BFS
        while(improved == 1) {
            if(threadIdx.x == 0) improved = 0;
            for(int node = threadIdx.x; node < n; node += blockDim.x) {
                for(int edge = rp[node]; edge < rp[node + 1]; edge++) {
                    int adj = ci[edge];
                    if(dist[(s * n) + adj] == level && dist[(s * n) + node] == -1) {
                        dist[(s * n) + node] = level + 1;
                        improved = 1;
                    }
                    if(dist[(s * n) + adj] == level && dist[(s * n) + node] == level + 1) {
                        sigma[(s * n) + node] += (float) sigma[(s * n) + adj];
                    }
                }
            }
            if(threadIdx.x == 0) level++;
            __syncthreads();
        }

        int dist_sum = 0;
        int dist2_cnt = 0;

        // DISTANCE ADDER
        if(threadIdx.x == 0) {
            for(int i = 0; i < n; i++) {
                if(dist[(s * n) + i] > 0) {
                    if(dist[(s * n) + i] <= 2) dist2_cnt++;
                    dist_sum += dist[(s * n) + i];
                }
            }
            results[n + s] = dist2_cnt; // degree 2
            results[2 * n + s] = (float) n / dist_sum; // closeness cent.
        }

	    while(level > 0) {
	        for(int node = threadIdx.x; node < n; node += blockDim.x) {
                if(dist[s * n + node] == level){
                    for(int edge = rp[node]; edge < rp[node + 1]; edge++) {
                        int adj = ci[edge];
                        if(dist[(s * n) + adj] + 1 == dist[(s * n) + node]) {
                            atomicAdd(&delta[(s * n) + adj], (sigma[(s * n) + adj] * 1.0) / sigma[(s * n) + node] * (1 + delta[(s * n) + node]));
                        }
                    }
                    atomicAdd(&results[3 * n + node], delta[(s * n) + node] / 2);
                }
            }
            if(threadIdx.x == 0) level--;
            __syncthreads();
	    }
    }
}

float* compute_centralities(int *rp, int *ci, int n, float &time_taken) {
    const int BLOCK_COUNT = MAX_MEMORY / (4 * 3 * n);
    int *sigma, *dist;
    float *delta, *d_results;

    cudaMalloc((void **)&d_results, sizeof(float) * n * 4);
    cudaMalloc((void **)&sigma, sizeof(int) * n * BLOCK_COUNT);
    cudaMalloc((void **)&dist, sizeof(int) * n * BLOCK_COUNT);
    cudaMalloc((void **)&delta, sizeof(float) * n * BLOCK_COUNT);

    cudaMemset(dist, -1, sizeof(int) * n * BLOCK_COUNT);
    cudaMemset(sigma, 0, sizeof(int) * n * BLOCK_COUNT);
    cudaMemset(delta, 0, sizeof(float) * n * BLOCK_COUNT);
    cudaMemset(d_results, 0, sizeof(float) * 4 * n);

    cudaEvent_t start, end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start);

    cent_kernel<<<BLOCK_COUNT, THREAD_COUNT>>>(d_results, dist, sigma, delta, rp, ci, n);
    cudaDeviceSynchronize();
    
    cudaEventRecord(end);
    cudaEventSynchronize(end);

    cudaEventElapsedTime(&time_taken, start, end);

    float *results;
    cudaMallocHost((void **)&results, sizeof(float) * n * 4);
    cudaMemcpy(results, d_results, sizeof(float) * n * 4, cudaMemcpyDeviceToHost);

    cudaFree(sigma);
    cudaFree(dist);
    cudaFree(delta);
    cudaFree(d_results);
    
    cudaDeviceSynchronize();
    return results;
}

int main()
{
    cudaSetDevice(0);

    std::string filename = "../data/wing_nodal.mtx";
    int *row_ptr, *col_ind;
    int num_nodes, num_edges;
    read_graph(filename, row_ptr, col_ind, num_nodes, num_edges);
    printf("[INFO] Graph is read: %s.\n", filename.c_str());

    int *rp;
    int *ci;

    cudaMalloc((void **)&rp, sizeof(int) * (num_nodes + 1));
    cudaMalloc((void **)&ci, sizeof(int) * (num_edges * 2));
    printf("[INFO] CUDA memory parameters are allocated for kernel function.\n");

    cudaMemcpy(rp, row_ptr, sizeof(int) * (num_nodes + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(ci, col_ind, sizeof(int) * (num_edges * 2), cudaMemcpyHostToDevice);

    printf("[INFO] CUDA memory parameters are set for kernel function.\n");
  
    float time_taken;  
    float *results = compute_centralities(rp, ci, num_nodes, time_taken);

    printf("[INFO] Kernel function is finished.\n");

    printf("Centrality Results:\n");
    for (int i = 0; i < num_nodes; i++)
    {
        printf("%.5f; %.5f; %.5f; %.5f\n ", results[i], results[num_nodes + i], results[2 * num_nodes + i], results[3 * num_nodes + i]);
    }
    printf("[INFO] Kernel call is ended in: %.5f ms.\n", time_taken);

    cudaFreeHost(results);
    cudaFreeHost(row_ptr);
    cudaFreeHost(col_ind);
    cudaFree(rp);
    cudaFree(ci);

    return 0;
}


