g++ -c BCTool.cpp
g++ -c Cell.cpp
g++ -c System.cpp
g++ -c CStemCell.cpp
g++ -c Random.cpp
g++ BCTool.o System.o CStemCell.o Cell.o Random.o -o bct_StemCell
