all:
	g++ sequence_container.cpp sequence.cpp main.cpp -std=c++11 -Wall -Wextra -pthread -lz -g -O3 -o pbclip
