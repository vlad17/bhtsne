all: sptree.cpp main.cpp tsne.hpp tsne.tpp
	g++ -std=c++11 sptree.cpp main.cpp -o bh_tsne -O3
