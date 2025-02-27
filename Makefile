CXX ?= g++
TARGET_EXEC ?= parsim

# paths
SRC_PATH := src
BUILD_PATH := build
BIN_PATH := $(BUILD_PATH)/bin

# sources
SOURCES := $(shell find $(SRC_PATH) -name *.cpp)
#OBJS := $(SOURCES:%=$(BIN_PATH)/%.o)
OBJS := $(patsubst %, $(BIN_PATH)/%.o, $(SOURCES))

# includes
INCLUDE_PATH := include
INCLUDE_FLAG := $(addprefix -I, $(INCLUDE_PATH))

# flags
CPPFLAGS ?= $(INCLUDE_FLAG) -O3
LDFLAGS ?= -fopenmp

# final target
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
# Default rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(BUILD_PATH)/%.o: $(SRC_PATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) -r $(BUILD_PATH)
	$(RM) $(TARGET_EXEC)
