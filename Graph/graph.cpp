#include "./graph.h"

int read_graphs(string &fname, int &num_nodes, int &num_edges, vector<int> &row_ptr, vector<int> &col_ind)
{
	ifstream input(fname.c_str());
	if (input.fail())
	{
		return -1;
	}
	// read graph
	string line = "%";
	while (line.find("%") != string::npos)
	{
		getline(input, line);
	}
	istringstream ss(line);
	ss >> num_nodes >> num_nodes >> num_edges;
	int v1, v2;
	double weight;

	vector<int> renameArr(num_nodes, -1);
	int counter = 0;
	bool eliminateUnused = true;

	vector<vector<int>> adj_list(num_nodes);
	for (int i = 0; i < num_edges; i++)
	{
		getline(input, line);
		istringstream inp(line);
		inp >> v1 >> v2;
		v1--; // make it 0 based
		v2--;

		//for detecting vetices that are unused
		if (renameArr[v1] == -1 && eliminateUnused)
		{
			renameArr[v1] = counter;
			v1 = counter;
			counter++;
		}
		else if (eliminateUnused)
		{
			v1 = renameArr[v1];
		}
		if (renameArr[v2] == -1 && eliminateUnused)
		{
			renameArr[v2] = counter;
			v2 = counter;
			counter++;
		}
		else if (eliminateUnused)
		{
			v2 = renameArr[v2];
		}

		if (v1 != v2)
		{
			adj_list[v1].push_back(v2); // add the edge v1->v2
			adj_list[v2].push_back(v1); // add the edge v2->v1
		}
	}
	if (eliminateUnused)
	{
		num_nodes = counter;
	}

	row_ptr = vector<int>(num_nodes + 1);
	col_ind = vector<int>(2 * num_edges);
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
	//cout << "nof nodes " << num_nodes << endl;
	//cout << "nof edges " << num_edges << endl;
	return 0;
}

void get_filenames(vector<string> &filenames, const vector<string> &locations)
{
    for (int i = 0; i < locations.size(); i++)
    {
        if (auto dir = opendir(locations[i].c_str()))
        {
            while (auto f = readdir(dir))
            {
                if (!f->d_name || f->d_name[0] == '.')
                    continue;

                string path = locations[i] + f->d_name;
                filenames.push_back(path);
            }
        }
    }
}

void bfs(int start_node, int num_nodes, vector<int> row_ptr, vector<int> col_ind, vector<int> &distance_arr, int step_size)
{
	vector<int> frontier(num_nodes, -1);
	frontier[0] = start_node;
	int queuestart = 0, queueend = 1, frontsize = 0;
	distance_arr.assign(num_nodes, -1); // every node is unvisited
	distance_arr[start_node] = 0;		// distance from a node to itself is 0
	int dist = 1;						// initial distance
	bool improvement = true;
	while (improvement && dist <= step_size)
	{
		improvement = false;
		do
		{
			int front = frontier[queuestart++];
			for (int edge = row_ptr[front]; edge < row_ptr[front + 1]; edge++)
			{ // for each adjacent of front
				int adj = col_ind[edge];
				if (distance_arr[adj] == -1)
				{ // if it is not visited
					improvement = true;
					frontier[queueend + frontsize++] = adj; // place it into next location (new frontier)
					distance_arr[adj] = dist;				// assign corresponding distance
				}
			}
		} while (queuestart < queueend);
		queueend += frontsize; // add the offset
		frontsize = 0;		   // reset the offset
		dist++;				   // next frontier will be further
	}
}

// for large scaled graphs, an approximating function of closeness centrality
void closeness_centrality_approx(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering, int size)
{
	if (size > num_nodes)
	{
		size = num_nodes;
	}
	vector<vector<int>> dist_arr(size, vector<int>(num_nodes));
	//	omp_set_num_threads(32);
	#pragma omp parallel for num_threads(32) schedule(dynamic)
	for (int v = 0; v < size; v++)
	{
		//srand(1);
		int start_node = rand() % num_nodes;
		bfs(start_node, num_nodes, row_ptr, col_ind, dist_arr[v]); // take distance array for node v
	}
	#pragma omp barrier

	#pragma omp parallel for num_threads(32)
	for (int v = 0; v < num_nodes; v++)
	{
		int sum_of_dist = 0;
		for (int i = 0; i < size; i++)
		{
			sum_of_dist += dist_arr[i][v];
		}
		float coeff = sum_of_dist > 0 ? (float)size / sum_of_dist : 0; // if coefficient is negative(meaning that graph is not connected) assign to 0
		ordering[v] = make_pair(v, coeff);
	}
}

void closeness_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering)
{
	vector<int> dist_arr;
	for (int v = 0; v < num_nodes; v++)
	{
		bfs(v, num_nodes, row_ptr, col_ind, dist_arr);							// take distance array for node v
		int sum_of_dist = std::accumulate(dist_arr.begin(), dist_arr.end(), 0); // sum of d(v,x) for all x in the graph
		float coeff = sum_of_dist > 0 ? (float)num_nodes / sum_of_dist : 0;		// if coefficient is negative(meaning that graph is not connected) assign to 0
		ordering[v] = make_pair(v, coeff);
	}
}

void degree_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering)
{
	//	omp_set_num_threads(32);
	#pragma omp parallel for num_threads(32) schedule(dynamic)
	for (int v = 0; v < num_nodes; v++)
	{
		ordering[v] = make_pair(v, row_ptr[v + 1] - row_ptr[v]);
	}
}

void degree2_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering)
{
	//	omp_set_num_threads(32);
	#pragma omp parallel for num_threads(8) schedule(dynamic)
	for (int v = 0; v < num_nodes; v++)
	{
		vector<int> dist_arr(num_nodes);
		bfs(v, num_nodes, row_ptr, col_ind, dist_arr, 2); // take distance array for node v
		int count = 0, val;
		for (int i = 0; i < dist_arr.size(); i++)
		{
			val = dist_arr[i];
			if (val == 1 || val == 2)
			{
				count++;
			}
		}
		ordering[v] = make_pair(v, count);
	}
}
