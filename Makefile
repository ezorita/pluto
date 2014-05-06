SRC_DIR= src
BIN_DIR= .
OBJ_DIR= obj

INDEX_SOURCES= plutoindex.c
INDEX_HEADERS= plutoindex.h plutocore.h
INDEX_OBJECTS= plutocore.o

QUERY_SOURCES= plutoquery.c
QUERY_HEADERS= plutocore.h
QUERY_OBJECTS= plutocore.o

INDEX_OBJ= $(addprefix $(OBJ_DIR)/,$(INDEX_OBJECTS))
INDEX_DEPS= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES) $(INDEX_HEADERS)) $(INDEX_OBJ)
INDEX_SRC= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES))

QUERY_OBJ= $(addprefix $(OBJ_DIR)/,$(QUERY_OBJECTS))
QUERY_DEPS= $(addprefix $(SRC_DIR)/,$(QUERY_SOURCES) $(QUERY_HEADERS)) $(QUERY_OBJ)
QUERY_SRC= $(addprefix $(SRC_DIR)/,$(QUERY_SOURCES))

INCLUDES= $(addprefix -I,$(SRC_DIR))
CFLAGS= -std=c99 -Wall -g -O3
LDLIBS= -lm
CC= gcc

all: plutoindex plutoquery

plutoindex: $(INDEX_DEPS)
	$(CC) $(CFLAGS) $(INDEX_SRC) $(INDEX_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

plutoquery: $(QUERY_DEPS)
	$(CC) $(CFLAGS) $(QUERY_SRC) $(QUERY_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(wildcard $(OBJ_DIR)/*) $(BIN_DIR)/plutoindex
