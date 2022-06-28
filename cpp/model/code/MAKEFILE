# a build command to build myprogram executable from myprogram.o and mylib.lib
all:main.cpp
	g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 main.cpp -Ofast -march=native -s -lgsl -lgslcblas
debug : main.cpp
	g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 main.cpp -g -lgsl -lgslcblas 
clean:
	$(RM) a.out
