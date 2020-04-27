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
	bool eliminateUnused = false;

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

void print_centralities(string filename, vector<vector<float>> &result) {
	cout << endl
		 << "-----------------------------------------------------" << endl
		 << "| Results for " << filename << endl
		 << "-----------------------------------------------------" << endl
		 << "|    Deg1    |    Deg2    |    Clos    |    Betw    |" << endl
		 << "    ------       ------       ------       ------   " << endl;

	for(int i = 0; i < result[0].size(); i++) {
		cout << "|";
		for(int j = 0; j < 4; j++) {
			string res = to_string(result[j][i]);
			if(result[j][i] >= 100) cout << " " <<  res << " |";
			else if(result[j][i] >= 10) cout << " " <<  res << "  |";
			else if(result[j][i] < 0) cout << " " <<  res << "  |";
			else cout << "  " <<  res << "  |";
		}
		cout << endl;
	}
}

/**
 * num_nodes			: number of nodes in the graph
 * row_ptr and col_ind	: constant graph data holders
 * result				: 1st dimension is the size 4 (number of closeness measures included), 2nd dimension is the 
 * 			 			  values of all nodes. Shape is = (4 x num_nodes). If one of the measures is not wanted, then 
 * 						  the values of all nodes for this measurement will be 0.
 * to_calculate			: boolean vector of size 4. 1st index for degree_centrality, 2nd for degree_2_centrality
 * 						  3rd for closeness_centrality, 4th for betweennes centrality
 * step_size			: for large-sized graphs an approximation parameter that indicates how many random nodes should
 * 						  be calculated to approximate (for closeness and betweenness centralities)
 * 						  !!! NOT IMPLEMENTED YET !!!
 */
void compute_centralities(const int &num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, 
						  vector<vector<float>> &result, const vector<bool> &to_calculate, int step_size) {
	
	// initial input check
	int num_centralities = 0;
	if(to_calculate.size() != 4) {
		printf("Error: correct number of calculation flags are not set!\n");
		return;
	}
	for(int i = 0; i < 4; i++) {
		if(to_calculate[i]) num_centralities++;
	}
	if(num_centralities == 0) {
		printf("Error: none of the centrality calculations are set true!\n");
		return;
	}
	// parameter initialization
	result = vector<vector<float>>(4, vector<float>(num_nodes, 0));

	#pragma omp parallel for shared(result, row_ptr, col_ind, to_calculate)
	for(int s = 0; s < num_nodes; s++) {
		vector<int> queue; 						// non-increasing order of nodes to be read will be hold
		queue.push_back(s);
		int front = 0;							// the index of last processed element will be hold
		vector<int> dist(num_nodes, -1); 		// distances of all nodes to node s
		int dist2_counter = 0;					// counts nodes with distances <= 2
		vector<vector<int>> pred(num_nodes);	// for each node, its predecessors are held
		vector<int> sigma(num_nodes, 0);		// number of shortest paths from s to the index
		sigma[s] = 1;
		
		// for degree 1 centrality, value calculation
		if(to_calculate[0]) result[0][s] = row_ptr[s + 1] - row_ptr[s];
		// if check for eliminating extra calculation if only degree 1 is requested
		if(!to_calculate[1] && !to_calculate[2] && !to_calculate[3]) continue;

		while(front < queue.size()) {  // BFS Phase
			int v = queue[front];
			front++;
			// if check for eliminating unnecessary computation if closeness and betweennes is not requested
			if(dist[v] == 2 && !to_calculate[2] && !to_calculate[3]) break;

			for(int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
				int w = col_ind[edge];
				// to calculate the distance of w to v
				if(dist[w] < 0) {
					dist[w] = dist[v] + 1;
					if(dist[w] <= 2) dist2_counter++;
					queue.push_back(w);
				}
				// to see if edge (v, w) in the shortest path
				if(to_calculate[3] && dist[w] == dist[v] + 1) {
					sigma[w] += sigma[v];
					pred[w].push_back(v);
				}
			}
		}
		// for degree 2 centrality, value calculation 
		if(to_calculate[1]) result[1][s] = dist2_counter;
		// for closeness centrality, value calculation 
		#pragma omp task shared(result, queue, dist)
		{
			if(to_calculate[2]) {
				int sum_of_dist = 0;
				for(int &item: dist){
					if(item != -1) sum_of_dist += item;
				}
				#pragma omp critical 
				{
					result[2][s] = (float) 1 / sum_of_dist;
				}
			}
		}

		// if check for eliminating extra calculation if betweenness is not requested
		#pragma omp task shared(result, queue, sigma, pred)
		{
			vector<float> delta(num_nodes, 0);		// dependency of s to index node
			if(to_calculate[3]) {
				for(int i = 0; i < queue.size(); i++) {
					int w = queue[i];
					for(int &v: pred[w]) {
						# pragma omp critical
						{
							delta[v] += (float) ((float) sigma[v] / sigma[w]) * (1 + delta[w]);
							if(w != s) result[3][w] += delta[w];
						}
					}
				}
			}
		}
		#pragma omp taskwait
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

void betweenness_centrality(int num_nodes, const vector<int> &row_ptr, const vector<int> &col_ind, vector<pair<int, float>> &ordering, int size)
{
	queue<int> Q;
	stack<int> S;
	ordering = vector<pair<int, float>> (num_nodes, {0, 0});
	for(int s = 0; s < num_nodes; s++) {
		vector<int> dist(num_nodes, -1);
		vector<vector<int>> pred(num_nodes);
		vector<int> sigma(num_nodes, 0);
		dist[s] = 0;
		sigma[s] = 1;
		Q.push(s);

		while(!Q.empty()) {
			int v = Q.front();
			Q.pop();
			S.push(v);
			for(int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
				const int &w = col_ind[edge];
				if(dist[w] < 0) {
					dist[w] = dist[v] + 1;
					Q.push(w);
				}
				if(dist[w] == dist[v] + 1) {		
					sigma[w] = sigma[w] + sigma[v];
					pred[w].push_back(v);
				}
			}
		}

		vector<float> delta(num_nodes, 0);
		while(!S.empty()) {
			int w = S.top();
			S.pop();
			for(int &v: pred[w]) {
				cout << sigma[v] << " - " << sigma[w] << " - " << delta[w] << endl;
				delta[v] += ((float) (sigma[v] / sigma[w])) * (1 + delta[w]);
			}
			if(w != s) {
				ordering[w].first = w;
				ordering[w].second += delta[w];
			}
		}
		cout << "In the end!" << endl;
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


