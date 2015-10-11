all: sptree.cpp main.cpp tsne.hpp tsne.tpp gdelt.hpp
	g++ -std=c++11 sptree.cpp main.cpp -o bh_tsne -O3
