.SILENT: run
run:
	g++ ./Graph/graph.h ./Graph/graph.cpp main.cpp -O3 -fopenmp
	./a.out 1
	rm ./a.out ./Graph/graph.h.gch