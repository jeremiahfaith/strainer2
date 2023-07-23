#ifndef BIO_HASH_H
#define BIO_HASH_H 1

/**
    \file BIO_hash.h
    \brief string hash library

  <pre><b>Example:</b>
        // initialize hash of size 100
        BIO_hash h=BIO_initHash(100);
	int x = 25;
	int y = 15;
	int *temp;

	// insert some stuff; format is ( <i>hash</i>, <i>key {must be a string}</i>, <i>pointer to data</i> )
        BIO_addHashData(h, "big", &x);
        BIO_addHashData(h, "small", &y);

	// you must cast back as it returns a pointer to void
        temp=(int *)BIO_searchHash(h, "big");
        printf("big is \%d", *temp);

        // returns a list of all keys
        keys=BIO_getHashKeys(h);
        for (i=0; i<BIO_getHashSize(h); i++)
                printf("key is \%s", keys[i]);

	// free space allocated to hash
        BIO_destroyHash(h);



  </pre>
*/



#ifdef __cplusplus
extern "C" {
#endif

typedef struct item {
        void *DATA;
        char *key;
} BIO_hashItem;


typedef struct hash {
       unsigned int M; // hash size
       unsigned int N; // number of items in hash
       BIO_hashItem *data;  // array of the hashed data
} *BIO_hash;


/** returns a hash. if \a initialSize = 0; 
 * it starts at capacity  DEFAULT_HASH_SIZE */
BIO_hash BIO_initHash(int initialSize);

/** inserts \a data into \a h, which can be retrieved using \a key */
void BIO_addHashData(BIO_hash h, const char *key, void *data);

/** inserts \a key into \a h, which can be retrieved using \a key; for use when only presence/absence of key is needed;
*   when BIO_searchHash() is called on this type of data it returns the key string */
void BIO_addHashKeyOnly(BIO_hash h, const char *key);

/** removes \a key from \a h. the data in the hash must be freed separately and AFTER calling remove key */
void* BIO_removeHashKey(BIO_hash h, const char *key);

/** returns a void pointer to the data in \a h corresponding to \a key;
 * returns NULL when the key does not exist in the hash */
void* BIO_searchHash(const BIO_hash h, const char *key);

/** returns an a pointer to a pointer of strings containing all
 * the keys; this array of character arrays must be freed
 * and can be done with BIO_destroyHashKeys(char **);
 * note that if you call BIO_hashAllKeys, then insert a new
 * item, then call BIO_getHashSize(\a h) to determine the proper number of keys
 * to cycle through, you will overstep the bounds of the array and crash your 
 * program; Note that after you call BIO_destroyHash() the keys are destroyed 
 * and the char** array returned by this function is no longer valid
 */
char** BIO_getHashKeys(const BIO_hash h); 
/** frees the keys allocated by BIO_getHashKeys() */
void BIO_destroyHashKeys(char **keys); 


/** sets the default hash size */
void BIO_setDefaultHashSize(int size);


/** frees all space used by \a h.
 * NOTE: it is the USERS responsibility to free any dynamically allocated data
 * placed in the hash     */
void BIO_destroyHash(BIO_hash h);
/** same as hash BIO_destroyHash(), but it frees any data in the hash as well
 *  (note this will ONLY work if you inserted a simple data structure that can be
 *  free with one call to free() )*/
void BIO_destroyHashD(BIO_hash h);

/** same as hash BIO_destroyHash(), but it frees any data in the hash as well
 *  by allowing the developer to pass a pointer to a function that will be called 
 *  to free the data for that value; this if more efficient than getting all the 
 *  keys by BIO_getHashKeys() and freeing the data by retrieving them one-at-a-time */ 
void BIO_destroyHashF(BIO_hash h, void *freeDATA(void *));

/** current size of BIO_hash \a h */
#define BIO_getHashSize(h) (h->N)
/** current capacity of BIO_hash \a h */
#define BIO_getHashCapacity(h) (h->M)


#ifdef __cplusplus
}
#endif

#endif // BIO_HASH_H
