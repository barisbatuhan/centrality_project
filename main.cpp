#include "./Graph/graph.h"

int main(int argc, char **argv)
{
    vector<string> locations = {"./data/small/"}; 
    vector<string> files;
    get_filenames(files, locations);
    //omp_set_num_threads(atoi(argv[1]));

    double total_time = 0;
    for (int i = files.size() - 1; i < files.size(); i++)
    {
        vector<int> row_ptr, col_ind;
        int num_nodes, num_edges;
        read_graphs(files[i], num_nodes, num_edges, row_ptr, col_ind);

        vector<vector<float>> centralities;
        vector<bool> to_calculate = {true, true, true, true};
        double start = omp_get_wtime();
        compute_centralities(num_nodes, row_ptr, col_ind, centralities, to_calculate);
        total_time += omp_get_wtime() - start;
        print_centralities(files[i], centralities);
        // cout << "Completed: " << files[i] << endl;
    }
    cout << "Num threads: " << omp_get_max_threads() 
         << " - Total Time for " << files.size() << " graphs: " << total_time  << " sec. " << endl;

    return 0;
}
