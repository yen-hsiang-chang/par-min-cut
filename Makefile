BUILD_DIR ?= build

GBBS_SRC = gbbs/gbbs/io.cc gbbs/gbbs/graph_io.cc
MIN_CUT_SRC = main.cpp

SRCS := $(MIN_CUT_SRC) $(GBBS_SRC)
OBJS := $(SRCS:%=$(BUILD_DIR)/obj/%.o)

CFLAGS = -O3 -std=c++17 -mcx16 -latomic -march=native -fopenmp -DNDEBUG -DPARLAY_USE_STD_ALLOC -DPARLAY_OPENMP

LFLAGS = -Iparlaylib/include -Igbbs/

min-cut: $(OBJS)
	CC $(CFLAGS) $(LFLAGS) $^ -o $(BUILD_DIR)/$@ 

$(BUILD_DIR)/obj/%.o: %
	mkdir -p $(@D)
	CC $(CFLAGS) $(LFLAGS) -c $^ -o $@

clean:
	-rm -rf $(BUILD_DIR)

.PHONY: clean
