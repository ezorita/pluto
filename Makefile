SRC_DIR= src
BIN_DIR= .
OBJ_DIR= obj

INDEX_SOURCES= plutoindex.c
INDEX_HEADERS= plutoindex.h plutocore.h
INDEX_OBJECTS= plutocore.o

QUERY_SOURCES= plutoquery.c
QUERY_HEADERS= plutocore.h
QUERY_OBJECTS= plutocore.o

ALIGN_SOURCES= plutoaligner.c mergesort.c sma.c
ALIGN_HEADERS= plutoaligner.h plutocore.h sma.h mergesort.h
ALIGN_OBJECTS= plutocore.o

INDEX_OBJ= $(addprefix $(OBJ_DIR)/,$(INDEX_OBJECTS))
INDEX_DEPS= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES) $(INDEX_HEADERS)) $(INDEX_OBJ)
INDEX_SRC= $(addprefix $(SRC_DIR)/,$(INDEX_SOURCES))

QUERY_OBJ= $(addprefix $(OBJ_DIR)/,$(QUERY_OBJECTS))
QUERY_DEPS= $(addprefix $(SRC_DIR)/,$(QUERY_SOURCES) $(QUERY_HEADERS)) $(QUERY_OBJ)
QUERY_SRC= $(addprefix $(SRC_DIR)/,$(QUERY_SOURCES))

ALIGN_OBJ= $(addprefix $(OBJ_DIR)/,$(ALIGN_OBJECTS))
ALIGN_DEPS= $(addprefix $(SRC_DIR)/,$(ALIGN_SOURCES) $(ALIGN_HEADERS)) $(ALIGN_OBJ)
ALIGN_SRC= $(addprefix $(SRC_DIR)/,$(ALIGN_SOURCES))

INCLUDES= $(addprefix -I,$(SRC_DIR))
CFLAGS= -std=c99 -Wall -g -O3
LDLIBS= -lm -lpthread
CC= gcc

all: plutoindex plutoquery pluto

pluto: $(ALIGN_DEPS)
	$(CC) $(CFLAGS) $(ALIGN_SRC) $(ALIGN_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

plutoindex: $(INDEX_DEPS)
	$(CC) $(CFLAGS) $(INDEX_SRC) $(INDEX_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

plutoquery: $(QUERY_DEPS)
	$(CC) $(CFLAGS) $(QUERY_SRC) $(QUERY_OBJ) $(LDLIBS) $(INCLUDES) -o $(BIN_DIR)/$@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(wildcard $(OBJ_DIR)/*) $(BIN_DIR)/plutoindex
