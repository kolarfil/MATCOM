

#############################################
#THE SIMPLEST VERSION OF MAKEFILE

# link commands

# matice : Matrix.o Vector.o Solver.o main.o IterativeSolver.o JacobiOR.o GaussSeidelOR.o SteepestDescent.o ConjugGrad.o PrecSteepestDescent.o PrecConjugGrad.o
# 	g++ Matrix.o Vector.o Solver.o IterativeSolver.o JacobiOR.o GaussSeidelOR.o SteepestDescent.o ConjugGrad.o PrecSteepestDescent.o PrecConjugGrad.o main.o -o matice


# #compilation commands

# Matrix.o : Matrix.cpp
# 	g++ -c -g  Matrix.cpp -o Matrix.o

# Vector.o : Vector.cpp
# 	g++ -c -g Vector.cpp -o Vector.o

# Solver.o : Solver.cpp
# 	g++ -c -g Solver.cpp -o Solver.o

# IterativeSolver.o : IterativeSolver.cpp
# 	g++ -c -g IterativeSolver.cpp -o IterativeSolver.o

# JacobiOR.o : JacobiOR.cpp
# 	g++ -c -g JacobiOR.cpp -o JacobiOR.o

# GaussSeidelOR.o : GaussSeidelOR.cpp
# 	g++ -c -g GaussSeidelOR.cpp -o GaussSeidelOR.o

# SteepestDescent.o : SteepestDescent.cpp
# 	g++ -c -g SteepestDescent.cpp -o SteepestDescent.o

# ConjugGrad.o : ConjugGrad.cpp
# 	g++ -c -g ConjugGrad.cpp -o ConjugGrad.o

# PrecSteepestDescent.o : PrecSteepestDescent.cpp
# 	g++ -c -g PrecSteepestDescent.cpp -o PrecSteepestDescent.o

# PrecConjugGrad.o : PrecConjugGrad.cpp
# 	g++ -c -g PrecConjugGrad.cpp -o PrecConjugGrad.o

# main.o : main.cpp
# 	g++ -c -g main.cpp -o main.o



###############################################



###############################################
#MAKEFILE USING PATTERN RULES


# CXX=g++
# CXXFLAGS=-g
# SOURCES=Matrix.cpp Vector.cpp Solver.cpp IterativeSolver.cpp JacobiOR.cpp GaussSeidelOR.cpp SteepestDescent.cpp ConjugGrad.cpp PrecSteepestDescent.cpp PrecConjugGrad.cpp main.cpp
# OBJECTS=$(SOURCES:.cpp=.o)


# matice: $(OBJECTS)
# 	$(CXX) $(inputs) -o $(output)

# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) -c $(input) -o $(output)


################################################

#Diferent try


CXX=g++
CXXFLAGS=-c -g
OBJDIR=../obj
BUILDDIR=../build
SOURCES=Matrix.cpp Vector.cpp Solver.cpp IterativeSolver.cpp JacobiOR.cpp GaussSeidelOR.cpp SteepestDescent.cpp ConjugGrad.cpp PrecSteepestDescent.cpp PrecConjugGrad.cpp main.cpp
_OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_OBJECTS))
EXECUTABLE=matice

all: $(SOURCES) $(BUILDDIR)/$(EXECUTABLE)
#all:$(BUILDDIR)/

#.cpp.o:
$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ 

$(BUILDDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ 


