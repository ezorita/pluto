SRC_DIR= src
BIN_DIR= .
OBJ_DIR= obj

INDEX_SOURCES= plutoindex.c
INDEX_HEADERS= plutoindex.h
INDEX_OBJECTS= 

INDEX_DEPS= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES) $(INDEX_HEADERS) $(INDEX_OBJECTS))
INDEX_SRC= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES))
INDEX_OBJ= $(addprefix $(OBJ_DIR)/,$(INDEX_OBJECTS))

INCLUDES= $(addprefix -I,$(SRC_DIR))
CFLAGS= -std=c99 -Wall -g -O3
LDLIBS= -lm
CC= gcc

all: plutoindex

plutoindex: $(INDEX_DEPS)
	$(CC) $(CFLAGS) $(INDEX_SRC) $(INDEX_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

clean:
	rm -f $(wildcard $(OBJ_DIR)/*) $(BIN_DIR)/plutoindex
