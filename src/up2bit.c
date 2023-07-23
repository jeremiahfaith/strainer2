#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "up2bit.h"

/* 
	TO DO

	1) switch to tommy hash (fixed size) to see if faster

*/

static const char TWO_BIT_TO_DNA[] = {'A', 'C', 'T', 'G'};

// http://stackoverflow.com/questions/111928/is-there-a-printf-converter-to-print-in-binary-format
void showbits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i=size-1;i>=0;i--)
    {
        for (j=7;j>=0;j--)
        {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}


/*
int main(int argc, char *argv[]) {
	char nucs[] = "CATGATACTACTGACTGCTATCGACTCTAC";
	uint64_t DNA = 0;
	char decode[33];

	for (int i=0; i<10; i++) {
		DNA = encode_DNA_2_bit(nucs, strlen(nucs));
		decode_DNA_2_bit(DNA, strlen(nucs), decode);
	}



	return 0;
}
*/


uint64_t encode_DNA_2_bit(const char *DNA, const int len) {
	int i;
	uint64_t up2bit = 0;
	char n;
//	printf("s:%s\n", DNA);
	for (i=0; i<len; i++) {
		up2bit <<= 2; // shift the DNA by 2 bits (only storing ACTG)
		n = DNA[i];
		n &= 0x6;
		n >>= 1;
		up2bit += n;
//		printf("%c\t%d\n", DNA[i], n);
//		showbits(sizeof(uint64_t), &up2bit);
	}

	return up2bit;
}


void decode_DNA_2_bit(uint64_t up2bit, const int len, char *decode) {
	int maxnucs = 32;
	int i;
	uint64_t front_mask = 0xC000000000000000;
	uint64_t n64;
	// remove the 0 bits at the beginning	
//	printf("mask\n");
//	showbits(sizeof(uint64_t), &front_mask);
//	printf("\ndecode\n");
///	showbits(sizeof(uint64_t), &up2bit);
	up2bit <<= ((maxnucs-len)*2);
//	printf("push bits to front\n");
//	showbits(sizeof(uint64_t), &up2bit);

	for (i=0; i<len; i++) {
		n64 = up2bit & front_mask;
		//showbits(sizeof(uint64_t), &n64);
		n64 >>= 62; // IS IT FASTER TO PULL THESE OFF IN REVERSE ORDER AND THEN REVERSE THE STRING?
//		printf("%d\n", n64);
		//showbits(sizeof(uint64_t), &n64);
//		showbits(sizeof(uint64_t), &up2bit);
		up2bit <<= 2;
		decode[i] = TWO_BIT_TO_DNA[n64]; // IS A SWITCH FASTER THAN A LOOKUP?
	}
	decode[len] = '\0';
//	printf("d:%s\n", decode);
}



