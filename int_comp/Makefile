run: build 
	./main < key.txt

build: 
	g++ crt.cpp fft_mult.cpp general.cpp ntru.cpp fheaes.cpp CtxtBit.cpp HE_Integer.cpp HE_Utils.cpp main.cpp -o main -lntl
clean:
	rm -rf main
