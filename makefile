.SILENT: run
run:
	g++ ./Graph/ugraph.h main.cpp -O3 -fopenmp
#	./a.exe 1
#	del .\a.exe .\Graph\ugraph.h.gch 	
	./a.out 1
	rm ./a.out ./Graph/ugraph.h.gch 