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

        vector<vector<pair<int, float>>> orders(3, vector<pair<int, float>>(num_nodes));
        degree_centrality(num_nodes, row_ptr, col_ind, orders[0]);
        degree2_centrality(num_nodes, row_ptr, col_ind, orders[1]);
        closeness_centrality(num_nodes, row_ptr, col_ind, orders[2]);

        cout << "Orders for " << files[i] << " is calculated..." << endl;
    }

    return 0;
}