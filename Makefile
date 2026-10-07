CC ?= cc
CFLAGS ?= -O3 -fPIC -std=c99 -Wall -Wextra
SRC_DIR := modules
LIB_DIR := lib
SRC := $(SRC_DIR)/spinos_kepler.c

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    SHARED_EXT := dylib
    SHARED_FLAGS := -dynamiclib
    LIB_NAME := libspinos_kepler.$(SHARED_EXT)
else ifeq ($(OS),Windows_NT)
    SHARED_EXT := dll
    SHARED_FLAGS := -shared
    LIB_NAME := spinos_kepler.$(SHARED_EXT)
else
    SHARED_EXT := so
    SHARED_FLAGS := -shared
    LIB_NAME := libspinos_kepler.$(SHARED_EXT)
endif

TARGET := $(LIB_DIR)/$(LIB_NAME)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	mkdir -p $(LIB_DIR)
	$(CC) $(CFLAGS) $(SHARED_FLAGS) -o $@ $< -lm

clean:
	rm -f $(LIB_DIR)/libspinos_kepler.dylib $(LIB_DIR)/libspinos_kepler.so $(LIB_DIR)/spinos_kepler.dll
