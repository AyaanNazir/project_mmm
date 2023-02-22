#If you have installed BLIS in a directory other than your home directory, set the HOME accordingly
#HOME      := /path_to_blis

-include objs.mk

# Make sure you have BLIS installed in your home directory
BLIS_LIB  := $(HOME)/blis/lib/libblis.a
BLIS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -I$(BLIS_INC) -m64 -mavx2 -std=c99 -march=native -fopenmp -D_POSIX_C_SOURCE=200809L
FFLAGS     := $(CFLAGS) 

PSIZE := -DNREPEATS=3   \
         -DP_BEGIN=48 \
         -DP_END=500 \
         -DP_INC=48

LDFLAGS    := -lpthread -m64 -lm -fopenmp

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o
SRC_OBJS   := driver.o gemm.o

# --------------------

.DEFAULT_TARGET := driver

debug: CFLAGS+=-DDEBUG 
debug: driver

driver: $(SRC_OBJS) $(MY_OBJS) $(UTIL_OBJS)
	$(LINKER) $(SRC_OBJS) $(MY_OBJS) $(UTIL_OBJS) $(BLIS_LIB) -o driver_gemm.x $(LDFLAGS) 

run: driver
	./driver_gemm.x  

%.o: %.c
	$(CC) $(CFLAGS) $(PSIZE) -c $< -o $@

# ---------------------                                                               
clean:
	rm -f *.o *~ core *.x
