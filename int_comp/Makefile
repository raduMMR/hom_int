run: build 
	./main < hNTRU/key.txt

build: 
	g++ -g hNTRU/crt.cpp hNTRU/fft_mult.cpp hNTRU/general.cpp hNTRU/ntru_multikey.cpp test_multikey.cpp main.cpp -o main -lntl
clean:
	rm -rf main