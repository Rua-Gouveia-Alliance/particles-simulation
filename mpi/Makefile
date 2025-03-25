CXX ?= mpicxx  # Use MPI compiler wrapper
TARGET_EXEC ?= parsim-mpi

# paths
SRC_PATH := src
BUILD_PATH := build
BIN_PATH := $(BUILD_PATH)/bin

# sources
SOURCES := $(shell find $(SRC_PATH) -name *.cpp)
OBJS := $(patsubst $(SRC_PATH)/%, $(BIN_PATH)/%.o, $(SOURCES))

# includes
INCLUDE_PATH := include
INCLUDE_FLAG := $(addprefix -I, $(INCLUDE_PATH))

# flags
CPPFLAGS ?= $(INCLUDE_FLAG) -O2
LDFLAGS ?= -fopenmp  # Keep OpenMP if needed, MPI is handled by the wrapper

# final target
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# build all dependencies
$(BIN_PATH)/%.cpp.o: $(SRC_PATH)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) -r $(BUILD_PATH)
	$(RM) $(TARGET_EXEC)

.PHONY: clean