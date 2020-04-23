#include "./Graph/graph.h"

int main()
{
    vector<string> locations = {"./data/small/"}; 
    vector<string> files;
    get_filenames(files, locations);

    for (int i = 0; i < files.size(); i++)
    {
        vector<int> row_ptr, col_ind;
        int num_nodes, num_edges;
        read_graphs(files[i], num_nodes, num_edges, row_ptr, col_ind);

        vector<vector<float>> centralities;
        vector<bool> to_calculate = {true, true, true, true};
        compute_centralities(num_nodes, row_ptr, col_ind, centralities, to_calculate);
        // print_centralities(files[i], centralities);
        cout << "Completed: " << files[i] << endl;
    }

    return 0;
}