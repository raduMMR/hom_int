run: build 
	./main < hNTRU/key.txt

build: 
	g++ hNTRU/crt.cpp hNTRU/fft_mult.cpp hNTRU/general.cpp hNTRU/ntru.cpp hNTRU/fheaes.cpp hcomp/CtxtBit.cpp hcomp/HE_Integer.cpp hcomp/HE_Utils.cpp test_functions.cpp main.cpp -o main -lntl
clean:
	rm -rf main
