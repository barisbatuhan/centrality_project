#include "./Graph/ugraph.h"
using namespace std;

int main(int argc, char **argv)
{
    vector<string> locations = {"./data/medium/"}; 
    // vector<string> files = {"./data/medium/epb1.mtx"};
    vector<string> files;
    UGraph::get_filenames(files, locations);
    omp_set_num_threads(atoi(argv[1]));

    double total_time = 0;
    for (unsigned int i = 0; i < files.size(); i++)
    {
        UGraph g(files[i]);
        vector<vector<float>> centralities;
        vector<bool> to_calculate = {true, true, true, true};
        double start = omp_get_wtime();
        g.compute_centralities(centralities, to_calculate);
        total_time += omp_get_wtime() - start;
        // UGraph::print_centralities(files[i], centralities);
    }
    cout << "Num threads: " << omp_get_max_threads() 
         << " - Total Time for " << files.size() << " graphs: " << total_time  << " sec. " << endl;

    return 0;
}