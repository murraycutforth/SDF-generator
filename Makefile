TARGET := testAPSS.exe
CFLAGS := -Wall -c -std=c++11 -g -O3
SOURCES := $(shell find . -type f -name '*.cpp')
OBJECTS := $(SOURCES:.cpp=.o)


$(TARGET) : $(OBJECTS)
	g++ $^ -O3 -o $(TARGET)

%.o: %.cpp
	g++ $(CFLAGS) -o $@ $<

clean:
	rm *.o $(TARGET)
