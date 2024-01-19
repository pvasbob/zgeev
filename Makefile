ALL = DiagMat.o test.o
CXX = g++
FLAG = -llapack -lblas

TARGET = test

$(TARGET): $(ALL)
	$(CXX) $(ALL) -o $(TARGET) $(FLAG)
