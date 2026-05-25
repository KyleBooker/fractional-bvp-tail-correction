# Build the main solver and the four tail-condition scheme variants.
# All binaries read homoclinic.in from the current directory.

CC      ?= gcc
CFLAGS  ?= -O2 -Wall
LDLIBS  ?= -lm

BINS = homoclinic \
       first_order_first_derivative \
       first_order_second_derivative \
       second_order_first_derivative \
       second_order_second_derivative

.PHONY: all clean

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BINS) $(addsuffix .exe,$(BINS)) homoclinic.out ic.in
