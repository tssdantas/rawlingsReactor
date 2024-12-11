all:
	g++-13 -std=c++14 -I/opt/homebrew/Cellar/boost/1.86.0_2/include -L/opt/homebrew/Cellar/boost/1.86.0_2/lib -lboost_system -Wall -g etano.cpp -o etano