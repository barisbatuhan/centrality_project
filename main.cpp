#include "./Graph/graph.h"

int main(int argc, char **argv)
{
    vector<string> locations = {"./data/large/"};
    vector<string> files;
    get_filenames(files, locations);
    //omp_set_num_threads(atoi(argv[1]));

    double total_time = 0;
    for (int i = files.size() - 1; i < files.size(); i++)
    {
        
        vector<vector<int>> adj_list;
        vector<int> row_ptr, col_ind;
        int num_nodes, num_edges;
        read_graphs(files[i], num_nodes, num_edges, row_ptr, col_ind);
        read_adjecancy(files[i], num_nodes, num_edges, adj_list);
        vector<vector<float>> centralities;
        vector<float> betweenCentrality(num_nodes, 0.0);
        vector<bool> to_calculate = {false, false, false, true};
        double start = omp_get_wtime();
        compute_centralities(num_nodes, row_ptr, col_ind, centralities, to_calculate);
        //betweennessCentrality(adj_list, betweenCentrality);
        total_time += omp_get_wtime() - start;
        print_centralities(files[i], centralities);
        /*
        for(int i=0; i<betweenCentrality.size(); i++){
            cout << betweenCentrality[i] << endl;
        }
         */
        // cout << "Completed: " << files[i] << endl;
    }
    cout << "Num threads: " << omp_get_max_threads() 
         << " - Total Time for " << files.size() << " graphs: " << total_time  << " sec. " << endl;

    return 0;
}
