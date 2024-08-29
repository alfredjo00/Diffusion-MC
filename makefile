CC = gcc

CFLAGS = \
	 -g \
	 -Iinclude \
	 -Wall \
	 -Werror \
	 -pedantic \
	 -fsanitize=address \
	 -fno-omit-frame-pointer

CFLAGS_OPT = \
	     -O3 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
      obj/run.o  \
      obj/dynamic_array.o \
      obj/diffusion_monte_carlo.o
	  
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

H3: obj _H3

_H3: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H3 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test
	rm -f H3

.PHONY: clean
