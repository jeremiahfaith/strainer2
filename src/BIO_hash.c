#include "BIO_hash.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MINIMUM_HASH_SIZE 10
#define CONST_B 27183

int DEFAULT_HASH_SIZE = 1000;

inline int hashU(const char *s, const int M);
void expand(BIO_hash);

BIO_hash BIO_initHash(int size) {
	//int i=0;
	BIO_hash h = (BIO_hash)malloc(sizeof(*h));

	if (!size)
		size=DEFAULT_HASH_SIZE;
	else if (size < MINIMUM_HASH_SIZE) 
		size=MINIMUM_HASH_SIZE;

	h->M = size;
	h->N=0;
	h->data = (BIO_hashItem*)calloc(sizeof(*h->data), h->M);

	/* this is slower but more platform independent 
	h->data = (BIO_hashItem*)malloc(sizeof(*h->data) * h->M);

	for (i=0; i<h->M; i++) {
		h->data[i].DATA=NULL;
		h->data[i].key=NULL;
	}
	*/
		
	return h;
}

void expand(BIO_hash h) {
	int i, M;
	BIO_hashItem *temp=h->data;
	M=h->M;
	h->M += h->M; // double storage capacity
	h->N = 0; // reset N, because we are going to refill the array
	h->data = (BIO_hashItem*)calloc(sizeof(*h->data), h->M);
	/* platform independent version, doesn't assume 0-bitwise == NULL
	h->data=(BIO_hashItem*)malloc(sizeof(BIO_hashItem) * h->M);
	for (i=0; i<h->M; i++) {
		h->data[i].DATA=NULL;
		h->data[i].key =NULL;
	}
	*/

	for (i=0; i<M; i++)
		if (temp[i].DATA != NULL) {
			BIO_addHashData(h, temp[i].key, temp[i].DATA);
			free(temp[i].key);
		}

	free(temp);	
}

void BIO_destroyHash(BIO_hash h) {
	int i;	
	/* free the local keys */
	for (i=0; i<BIO_getHashCapacity(h); i++) {
		if (h->data[i].DATA != NULL) 
			free(h->data[i].key);
	}

	/* free the BIO_hashItem array */
	free(h->data);
	/* free the hash */
	free(h);
}

void BIO_destroyHashD(BIO_hash h) {
	int i;	

	/* free the local keys */
	for (i=0; i<BIO_getHashCapacity(h); i++) {
		if (h->data[i].DATA != NULL) {
			free(h->data[i].key); 
			free(h->data[i].DATA);
		}
	}

	/* free the BIO_hashItem array */
	free(h->data);
	/* free the hash */
	free(h);
}

void BIO_destroyHashF(BIO_hash h, void *freeDATA(void *)) {
	int i;	

	/* free the local keys */
	for (i=0; i<BIO_getHashCapacity(h); i++) {
		if (h->data[i].DATA != NULL) {
			free(h->data[i].key); 
			freeDATA(h->data[i].DATA);
		}
	}

	/* free the BIO_hashItem array */
	free(h->data);
	/* free the hash */
	free(h);
}

void BIO_addHashKeyOnly(BIO_hash h, const char *key) {
	int i=hashU(key, h->M);
	char *data_and_key = strdup(key);
	while(h->data[i].DATA != NULL)
		i=(i+1) % h->M;

	h->data[i].key = data_and_key;
	h->data[i].DATA = data_and_key;
//	h->data[i].key = strdup(key);
//	h->data[i].DATA = data;

	/* make hash bigger if it is getting too big */
	if (h->N++ >= h->M/2) expand(h);
}




void BIO_addHashData(BIO_hash h, const char *key, void *data) {
	int i=hashU(key, h->M);
	while(h->data[i].DATA != NULL)
		i=(i+1) % h->M;

	h->data[i].key = strdup(key);
	h->data[i].DATA = data;

	/* make hash bigger if it is getting too big */
	if (h->N++ >= h->M/2) expand(h);
}

void* BIO_removeHashKey(BIO_hash h, const char *key) {
	int i=hashU(key, h->M);	
	void *data;
	while (h->data[i].DATA != NULL) {
		// see if keys match
		if (strcmp(key, h->data[i].key) == 0) {
			data = h->data[i].DATA; // grab the data to return one last time
			free(h->data[i].key);
			h->data[i].DATA = NULL;
			h->N--; // decrement the number of keys in hash
			return data;
		}
		else {
			i=(i+1) % h->M;
		}
	}

	return NULL;
}

void* BIO_searchHash(const BIO_hash h, const char *key) {
	int i=hashU(key, h->M);	
	while (h->data[i].DATA != NULL) {
		// see if keys match
		if (strcmp(key, h->data[i].key) == 0)
			return h->data[i].DATA;
		else
			i=(i+1) % h->M;
	}

	return NULL;
}

char** BIO_getHashKeys(const BIO_hash h) {
	int i;
//	char **keys = (char**)malloc(sizeof(char*) * (h->N+1));
	char **keys = (char**)malloc(sizeof(char*) * h->N);
	int count=0;
	for (i=0; i<h->M; i++)
		if (h->data[i].DATA != NULL)
			keys[count++]=h->data[i].key; /* don't store things multiple times to save memory */
	
	/* old; previous version copied the entire set of strings */
	// a little trick so I can deallocate the array 
//	keys[count]=NULL;

	return keys;
}

void BIO_destroyHashKeys(char **keys) {
//	int i=0;
//	while (keys[i]!=NULL)
//		free(keys[i++]);

	free(keys);
}


/*
 copied from glib/ghash.c
https://github.com/GNOME/glib/blob/master/glib/ghash.c
 * This function implements the widely used "djb" hash apparently
 * posted by Daniel Bernstein to comp.lang.c some time ago.  The 32
 * bit unsigned hash value starts at 5381 and for each byte 'c' in
 * the string, is updated: `hash = hash * 33 + c`. This function
 * uses the signed value of each byte.
*/
int hashU(const char *s, const int M) {
	const char *p;
	unsigned int h = 5381;

	for (p = s; *p != '\0'; p++)
		h = (h << 5) + h + *p;

	return h%M;
}

void BIO_setDefaultHashSize(int size) 
{ DEFAULT_HASH_SIZE = size; }


