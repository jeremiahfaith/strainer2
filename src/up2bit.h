
#ifndef UP2BIT_H
#define UP2BIT_H 1

/** 
 \file  up2bit.h
 \brief conversion of DNA (2-bit) 64-bit integer
 
 currently cannot handle more than 32 nucleotides (might add in the future using vector of uint64_t)
 ideas from => https://github.com/JohnLonginotto/ACGTrie/blob/master/docs/UP2BIT.md

 */

#ifdef __cplusplus
extern "C" {
#endif

uint64_t encode_DNA_2_bit(const char *DNA, const int len);
void decode_DNA_2_bit(uint64_t up2bit, const int len, char *nucs);



#ifdef __cplusplus
}
#endif

#endif // UP2BIT_H
