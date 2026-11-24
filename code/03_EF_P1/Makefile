CC      = gcc
CFLAGS  = -O3 -std=c99 -Wall -Wextra -Wpedantic
LDFLAGS = -lm

TARGET  = ef_p1_lshape
SRC     = ef_p1_lshape.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.dat *.pdf *.gif *.aux *.log *.out *.toc

.PHONY: all run clean
