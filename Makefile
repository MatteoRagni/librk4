CC = gcc
CFLAGS = -Wall -I. librk4.c

doc:
	doxygen .doxygen

test:
	$(CC) $(CFLAGS) test/test_ode_1.c -o test/test_ode_1
	$(CC) $(CFLAGS) test/test_ode_2.c -o test/test_ode_2
	$(CC) $(CFLAGS) test/test_ode_3.c -o test/test_ode_3

clean:
	rm -f test/test_ode_1 test/test_ode_2 test/test_ode_3
