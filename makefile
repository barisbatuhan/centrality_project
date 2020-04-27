.SILENT: run
run:
	g++ ./Graph/graph.h ./Graph/graph.cpp main.cpp -O3 -fopenmp
	./a.out
	rm ./a.out ./Graph/graph.h.gch
	
