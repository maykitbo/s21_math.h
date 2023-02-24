CC=gcc
SRC=s21_math.c
OBJ=s21_math.o
CFLAGS=-Wall -Werror -Wextra -std=c11
TESTFLAGS=-lcheck -lm -lpthread -lgcov

all: test gcov_report s21_math.a

gcov_report: test
	gcovr -r. --html --html-details -o tests/coverage_report.html

test: lib.a s21_tests.o $(OBJ)
	$(CC) $(FLAGS) lib.a s21_tests.o $(OBJ) $(TESTFLAGS) -o s21_test
	./s21_test

lib.a: $(OBJ)
	ar rcs lib.a $(OBJ)

s21_math.a: s21_math_verter.o
	ar rcs s21_math.a $(OBJ)

s21_tests.o: s21_tests.c
	$(CC) $(CFLAGS) --coverage -c $^

s21_math.o: $(SRC) 
	$(CC) $(CFLAGS) --coverage -c $^

s21_math_verter.o: $(SRC)
	$(CC) $(CFLAGS) -c $^

clean:
	rm -rf *.o *.gc* *.a
	rm -f s21_test
	rm -rf tests/*

