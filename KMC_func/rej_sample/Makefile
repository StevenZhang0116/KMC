CXX = g++
CXXFLAGS = -std=c++17 -O3 -w
TARGETS = rejection_sampling_1d rejection_sampling_2d

SRC_DIR = .
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

.PHONY: all clean

all: $(BUILD_DIR) $(TARGETS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS): %: $(OBJS)
	$(CXX) $(CXXFLAGS) $(filter $(SRC_DIR)/$(@F).cpp,$(SRCS)) -o $@ -lstdc++fs

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
