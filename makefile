all: main

main: trabalho_3_v2.o 
	g++ trabalho_3_v2.o -o main

CXXFLAGS = -O3 -march=native
trabalho_3_v2.o: trabalho_3_v2.cpp
	g++ $(CXXFLAGS) -c trabalho_3_v2.cpp

clean: 
	rm -f *.o main.exe
