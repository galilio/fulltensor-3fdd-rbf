CC      = gcc
CFLAGS  = -O3 -std=c99 -Wall -Wextra -Wpedantic
LDFLAGS = -lm

TARGET  = vf_tpfa_lshape
SRC     = vf_tpfa_lshape.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.dat *.pdf *.gif *.aux *.log *.out *.toc

.PHONY: all run clean
