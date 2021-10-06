main_BBM.out : main_BBM.o basic_particle.o 
	g++ -O3 -std=c++11 -fpermissive -o main_BBM.out basic_particle.o main_BBM.o -lgsl -lgslcblas -lm
main_BBM.o : main_BBM.cpp basic_particle.h
	g++ -O3 -std=c++11 -fpermissive -o main_BBM.o -c main_BBM.cpp
basic_particle.o : basic_particle.h
	g++ -O3 -c -Wall -std=c++11 -o basic_particle.o basic_particle.cpp -lgsl -lgslcblas -lm

