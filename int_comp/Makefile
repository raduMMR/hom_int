run: build 
	./main < Hom_NTRU_library/key.txt

build: 
	g++ Hom_NTRU_library/crt.cpp Hom_NTRU_library/fft_mult.cpp Hom_NTRU_library/general.cpp Hom_NTRU_library/ntru.cpp Hom_NTRU_library/fheaes.cpp CtxtBit.cpp HE_Integer.cpp HE_Utils.cpp main.cpp -o main -lntl
clean:
	rm -rf main
