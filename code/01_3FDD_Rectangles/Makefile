CC      = gcc
CFLAGS  = -O3 -std=c99 -Wall -Wextra -Wpedantic
LDFLAGS = -lm

TARGET  = tfdd_rect_lshape
SRC     = tfdd_rect_lshape.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.dat *.pdf *.aux *.log *.out *.toc

.PHONY: all run clean
