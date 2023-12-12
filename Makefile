# Makefile

# Compiler and flags
CXX = /usr/bin/gcc
CXXFLAGS = -Wall -Wextra -std=c++11 -arch arm64

# Source files
SRCS = main.cpp utils.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = my_program

# Build rule
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
