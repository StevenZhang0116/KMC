CXX = g++
CXXFLAGS = -std=c++17 -O3 -w
TARGETS = test_main

SRC_DIR = .
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
HDRS = $(wildcard $(SRC_DIR)/*.hpp)
HDRS += rejsample.hpp
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
OBJS += $(patsubst $(SRC_DIR)/%.hpp,$(BUILD_DIR)/%.o,$(HDRS))

.PHONY: all clean

all: $(BUILD_DIR) $(TARGETS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS): %: $(OBJS)
	$(CXX) $(CXXFLAGS) $(filter $(SRC_DIR)/$(@F).cpp,$(SRCS)) -o $@ -lstdc++fs

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
