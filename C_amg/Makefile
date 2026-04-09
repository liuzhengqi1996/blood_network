CC = gcc
CFLAGS = -O3 -Wall -Wextra -fPIC
LDFLAGS = -lm

SRCDIR  = src
INCDIR  = include
LIBDIR  = lib

SRCS = vector_ops.c sparse_ops.c io.c \
       build_hierarchy.c amg_cycle.c ichol.c amg_ichol_precond.c amg_pcg.c \
       rs_classical.c
OBJS = $(addprefix $(SRCDIR)/,$(SRCS:.c=.o))
HEADERS = $(INCDIR)/agg_lib.h $(INCDIR)/amg.h

PROG = test_amg
LIB  = $(LIBDIR)/libamg.so

all: $(PROG) libamg

$(PROG): test_amg.o libamg
	$(CC) $(CFLAGS) -o $@ $< -I$(INCDIR) -L$(LIBDIR) -lamg -Wl,-rpath,'$$ORIGIN/$(LIBDIR)' $(LDFLAGS)

libamg: $(SRCDIR)/amg_wrapper.o $(OBJS)
	$(CC) -shared -o $(LIB) $^ $(LDFLAGS)

test_amg.o: test_amg.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCDIR) -c -o $@ $<

$(SRCDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCDIR) -c -o $@ $<

clean:
	rm -f test_amg.o $(SRCDIR)/*.o $(PROG) $(LIB)

.PHONY: all clean libamg
