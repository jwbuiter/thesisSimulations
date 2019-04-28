CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

RFEA: RFEA.o
	$(CC) -o RFEA RFEA.o $(LDFLAGS)

RFEA.o: RFEA.cpp
	$(CC) -c -o RFEA.o RFEA.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o RFEA
