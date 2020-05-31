#include "ugraph.h"
using namespace std;

int main(int argc, char **argv)
{
    string file = "../data/wing_nodal.mtx";
    omp_set_num_threads(atoi(argv[1]));

    double total_time = 0;
    UGraph g(file);
    vector<vector<float>> centralities;
    vector<bool> to_calculate = {true, true, true, true};
    
    double start = omp_get_wtime();
    g.compute_centralities(centralities, to_calculate);
    total_time += omp_get_wtime() - start;
    
    UGraph::print_centralities(file, centralities);
    cout << "Total time: " << total_time << " seconds" << endl;

    return 0;
}
