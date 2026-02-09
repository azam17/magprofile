CC ?= gcc
CFLAGS = -Wall -Wextra -Werror -O2 -std=c11 -Ilib -Isrc
LDFLAGS = -lz -lm -lpthread

SRC_DIR = src
LIB_DIR = lib
TEST_DIR = test
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
MAIN_OBJ = $(BUILD_DIR)/main.o
LIB_OBJS = $(filter-out $(MAIN_OBJ),$(OBJS))

TARGET = halalseq

TEST_SRCS = $(wildcard $(TEST_DIR)/test_*.c)
TEST_BINS = $(patsubst $(TEST_DIR)/%.c,$(BUILD_DIR)/%,$(TEST_SRCS))

SDL2_CFLAGS = $(shell sdl2-config --cflags 2>/dev/null)
SDL2_LIBS   = $(shell sdl2-config --libs 2>/dev/null)

GUI_DIR = gui
GUI_SRCS = $(GUI_DIR)/main_gui.c $(GUI_DIR)/gui_analysis.c $(LIB_DIR)/tinyfiledialogs.c
GUI_OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(notdir $(GUI_SRCS)))

.PHONY: all clean test gui

all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $@

# Test targets
$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.c $(LIB_OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(LIB_OBJS) $(LDFLAGS) -o $@

test: $(TEST_BINS)
	@echo "=== Running HalalSeq Tests ==="
	@pass=0; fail=0; \
	for t in $(TEST_BINS); do \
		name=$$(basename $$t); \
		printf "  %-30s " "$$name"; \
		if $$t > /dev/null 2>&1; then \
			echo "PASS"; \
			pass=$$((pass + 1)); \
		else \
			echo "FAIL"; \
			fail=$$((fail + 1)); \
		fi; \
	done; \
	echo "=== $$pass passed, $$fail failed ===";\
	[ $$fail -eq 0 ]

# --- GUI targets ---
$(BUILD_DIR)/main_gui.o: $(GUI_DIR)/main_gui.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(SDL2_CFLAGS) -Igui -c $< -o $@

$(BUILD_DIR)/gui_analysis.o: $(GUI_DIR)/gui_analysis.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(SDL2_CFLAGS) -Igui -c $< -o $@

$(BUILD_DIR)/tinyfiledialogs.o: $(LIB_DIR)/tinyfiledialogs.c | $(BUILD_DIR)
	$(CC) -Wall -O2 -c $< -o $@

gui: $(LIB_OBJS) $(GUI_OBJS)
	mkdir -p HalalSeq.app/Contents/MacOS HalalSeq.app/Contents/Resources
	$(CC) $(LIB_OBJS) $(GUI_OBJS) $(LDFLAGS) $(SDL2_LIBS) \
	    -framework Cocoa -o HalalSeq.app/Contents/MacOS/HalalSeq
	cp $(GUI_DIR)/Info.plist HalalSeq.app/Contents/

clean:
	rm -rf $(BUILD_DIR) $(TARGET) HalalSeq.app
