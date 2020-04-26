.SILENT: run
run:
	g++ ./Graph/graph.h ./Graph/graph.cpp main.cpp -O3 -fopenmp
	./a.exe
	del .\a.exe .\Graph\graph.h.gch