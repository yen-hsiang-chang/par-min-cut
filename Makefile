BUILD_DIR ?= build

GBBS_SRC = gbbs/gbbs/io.cc
MIN_CUT_SRC = main.cpp

SRCS := $(GBBS_SRC) $(MIN_CUT_SRC)
OBJS := $(SRCS:%=$(BUILD_DIR)/obj/%.o)

JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -L$(JEMALLOCLD)
CFLAGS = -O3 -DNDEBUG -std=c++17 -mcx16 -latomic -march=native \
				 -fopenmp -DPARLAY_USE_STD_ALLOC # -DPARLAY_OPENMP 
LFLAGS = -Iparlaylib/include -Igbbs/ $(JEMALLOC)

min-cut: $(OBJS)
	CC $(CFLAGS) $(LFLAGS) $^ -o $(BUILD_DIR)/$@ 

$(BUILD_DIR)/obj/%.o: %
	mkdir -p $(@D)
	CC $(CFLAGS) $(LFLAGS) -c $^ -o $@

clean:
	-rm -rf $(BUILD_DIR)

.PHONY: clean
