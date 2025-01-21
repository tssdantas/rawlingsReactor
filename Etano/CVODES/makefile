CC = g++
CFLAGS = -I/usr/local/include
LDFLAGS = -L/usr/local/lib -lsundials_cvodes -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense

stiff_ode: stiff_ode.cpp
	$(CC) $(CFLAGS) stiff_ode.cpp -o stiff_ode $(LDFLAGS)

clean:
	rm -f stiff_ode