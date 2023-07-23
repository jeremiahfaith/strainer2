#ifndef BIO_SEQUENCE_H
#define BIO_SEQUENCE_H 1

#include <string.h>

/** 
 \file  BIO_sequence.h
 \brief manipulation and exploration DNA and Amino Acid sequences
 
 BIO_sequence.h contains basic functions for manipulation of biological 
 sequences.  The bioSeqsStruct structure is primary datatype for most of the functions 
 in the bioCPU library, but is most frequently used by it's typedef BIO_sequences.
 It stores the sequences, basic information about the sequence, provides buffers
 to allow more 'string-like' manipulation of sequences making 
 parser creation much simpler.  Most of the functions in this particular header
 pertain to rudimentary sequence access.  Slightly more biological
 functions include BIO_reverseComplement(), BIO_getTranslation(), BIO_getConsensus,
 and BIO_getGapSites().

 */

#ifdef __cplusplus
extern "C" {
#endif


#define SEQ_NUM_CHAR_COUNTS 256
#define BIO_SEQ_BUFFER_SIZE 0

/** \brief biological sequence types */
typedef enum bioSeqType { 
	BIO_NUCLEOTIDE,    /**< nucleotide sequence */
	BIO_AMINO_ACID,    /**< amino acid sequence */
	BIO_UNKNOWN_TYPE,  /**< sequence type could not be inferred */
	BIO_UNCHECKED_TYPE /**< sequence type has not been set */
} BIO_seqType;

/** \brief genetic codes / translation tables. 
 *  These tables were taken from the list of genetic codes used on NCBI's
 *  website
 */
typedef enum bioGeneticCode { 
	BIO_STANDARD_CODE,          /**< code used by most organism */
	BIO_VERTEBRATE_MT_CODE,     /**< specific to vertebrate mitochondrial genomes */
	BIO_INVERTEBRATE_MT_CODE,   /**< specific to invertebrate mitochondrial genomes */
	BIO_YEAST_MT_CODE,          /**< specific to yeast mitochondrial genomes. Note that CGA and CGC 
				      are absent from this code and are represented here by X */
	BIO_USER_DEFINED_CODE       /**< use the translation table set by BIO_setUserDefinedCode() */
} BIO_geneticCode;
/** \brief consensus sequence types   */
typedef enum bioConsensusType { 
/** for nucleotides creates a consensus representing each
    sequence by its IUPAC code for example <pre>
fish   ACACG
duck   AGAGG
dog    GTAGG
</pre> would yield consensus <i>RBASG</i>

for amino acids I've made up my own letter code with the six
remaining letters of the alphabet (if there is a standard code 
<b>please let me know</b>)
<pre>
B = AIVLFPM  = hydrophobic/nonpolar
J = STYHCNQW = polar
U = DE       = negatively charged
O = KRH      = positively charged
Z = DEKRH    = charged
X = everything else
</pre>
*/
	BIO_REDUNDANCY, 
	BIO_MAJORITY_RULE, /**< keeps the most frequent nucleotide/amino acid for a column */
	BIO_STRICT /**< only keeps a letter if it is completely conserved otherwise puts a gap */
} BIO_consensusType;

/** \brief constant to mark any column with a gap when using BIO_getGapSites() */
#define BIO_ANY_GAP 0.0
/** \brief constant to mark columns that contain ONLY gaps when using BIO_getGapSites() */
#define BIO_ALL_GAPPED 1.0

/** \brief constant for the number of amino acids in proteins */
#define BIO_NUM_AMINO_ACIDS 20
/** \brief constant for the number of postively charged amino acids */
extern int BIO_NUM_POSITIVE_CHARGE;
/** \brief constant for the number of negatively charged amino acids */
extern int BIO_NUM_NEGATIVE_CHARGE;
/** \brief constant for the number of charged amino acids.
 
   (this is of course just the union of BIO_NUM_POSITIVE_CHARGE and BIO_NUM_NEGATIVE_CHARGE)
 */
extern int BIO_NUM_CHARGED;
/** \brief constant for the number of hydrophobic/nonpolar amino acids */
extern int BIO_NUM_NONPOLAR;
/** \brief constant for the number of polar amino acids */
extern int BIO_NUM_POLAR;


#include <stdio.h>


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                        DATA STRUCTURES                               %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/** 
 * \brief stores individual sequences and relevant information about them
 
 * typically library users will use only a BIO_sequences typedef 
 * which contains a BIO_seqItem typedef
 * but they may also want to use the List structure directly for more low
 * level manipulation. If future compatability is a concern you shouldn't 
 * touch this structure, as you can obtain everything you need with 
 * accessor functions 
 */
typedef struct bioSeqItem {
	char *seq;		/**< the sequence */
	int *seqQual;           /**< the (optional) sequence quality scores */
	char *name;		/**< the sequence name */
	char *shortName;	/**< an abbreviated name, adjust length with 
				  BIO_setShortNameLength() */
	struct bioSeqItem *next;/**< allows linked-list of sequences */
	void *DATA;		/**< allows any arbitrary data to be associated
	                             with the record (not typically used) */
	int seqLength;		/**< length of sequence */
	int bufferSize;		/**< size of current buffer (for concatenating, mainly during file reading) */
} *BIO_seqItem;


/** \brief provides an abstraction of biological sequence data.
     
  Information in the structure is not typically manipulated directly, using accessor
  functions instead.  If you are not a <b><i>bioCPU</i></b> developer and/or don't want
  to learn about the guts of the library, I recommend
  looking at functions like BIO_getSeq() rIather than looking here.
  
  It may be used to maintain one or multiple sequences and multiple
  sequences may take the form of an alignment 
 */
typedef struct bioSeqsStruct {
	int numSeqs;		/**< number of sequences */
	int hasMatrix;		/**< has sequence as an array (mainly for dealing with alignments) */
	int shortNameLength;	/**< length to use when making short names */
	BIO_seqType seqType;	/**< type of sequence BIO_AMINO_ACID, BIO_NUCLEOTIDE */
	BIO_seqItem *matrix;	/**< matrix of seqs */
	BIO_seqItem seqs;	/**< linked list of seqs */
} *BIO_sequences;


extern char NUCLEOTIDES[];
extern char COMPLEMENT[];
extern int BIO_charCounts[SEQ_NUM_CHAR_COUNTS];
//extern int NUC_TO_INT[255];
extern int BLOSUM_50[][20];


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %           FUNCTIONS FOR ALLOCATING AND CREATING SEQUENCES            %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/** allocates a BIO_sequences structure and returns a pointer to the space */
BIO_sequences BIO_initSequences();

/** allocates a BIO_seqItem and returns a pointer to the space */
inline BIO_seqItem BIO_initSeqItem();

/** returns a copy of \a s */
BIO_sequences BIO_copySequences(const BIO_sequences s);

/** returns a copy of \a sI */
BIO_seqItem BIO_copySeqItem(const BIO_seqItem sI);

/** concatenates \a seq to sI allocating memory if necessary; \a strLength
 * is the length of \a seq.
 */ 
inline void BIO_catStr(BIO_seqItem sI, const char *seq, int strLength);

/** adds a BIO_seqItem \a sI
 *  It returns the newly added sequence as a BIO_seqItem, as during parsing
 *  this is convenient if you need to concatenate things again to the same sequence
 *  using BIO_seqCat(). Generally you don't need the returned BIO_seqItem, since it is
 *  also placed in \a s.
 */
BIO_seqItem BIO_addSeqItem(BIO_sequences s, BIO_seqItem sI);


/** allocates and creates a BIO_seqItem \a name for sequences \a seq.  Note that
    name <i>must NOT be NULL</i>.*/
BIO_seqItem BIO_createSeqItem(const char *name, const char *seq);


/** takes for example AGTCGAR and puts: AGTCGAG and AGTCGAA into actgSequence
 * actgSequence must be freed!  */ 
BIO_seqItem BIO_iupacToACGT(BIO_seqItem actgSequence);

/** gets BIO_seqItem corresponding to \a name and returns it's pointer. 
 * if BOOLEAN is true and \a name doesn't exist it creates it */
inline BIO_seqItem BIO_getSeqItem(BIO_sequences S, const char *name, int BOOLEAN);

/** gets BIO_seqItem at \a row in a matrix of sequences \a S; if \a S doesn't exist 
    it returns NULL. */
inline BIO_seqItem BIO_getSeqItemByIndex(BIO_sequences S, int row);


/** used to minimize the space used by S, only useful when BIO_seqCat
 * or other functions using buffers are utilized; it allocates only enough
 * space to store the sequence and frees the buffers */ 
void BIO_seqClean(BIO_sequences S);

/** frees up any and all space used by \a S */
void BIO_destroySequences(BIO_sequences S);

/** frees up any and all space used by \a sI */
inline void BIO_destroySeqItem(BIO_seqItem sI);




/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %              FUNCTIONS FOR manipulating sequence strings             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/** reverses string s in place; algorithm comes form the K&R book
 * note this will fail if someone has used char *myString="ACTACTACGA"
 * to initialize the array; must use malloc or char myString[]="ACTACTACGA"
 */
void BIO_reverseString(char[]);

/** reverse complements a sequence.
 * Note this will fail if someone has used char *myString="ACTACTACGA"
 * to initialize the array; must use malloc or char myString[]="ACTACTACGA".
 * Also, if you want to preserve the original sequence you should make a copy of
 * if or reverse complement again to get back to the original.
 */
void BIO_reverseComplement(char[]);

/** make a string upper case; This function useful when searching within strings or creating databases of sequences to make sure 
all the sequences are the same case.
 * Note this will fail if someone has used char *myString="ACTACTACGA"
 * to initialize the array; must use malloc or char myString[]="ACTACTACGA".
 */
void BIO_stringToUpper(char[]);

/** make a string lower case; This function useful when searching within strings or creating databases of sequences to make sure 
all the sequences are the same case.
 * Note this will fail if someone has used char *myString="ACTACTACGA"
 * to initialize the array; must use malloc or char myString[]="ACTACTACGA".
 */
void BIO_stringToLower(char[]);





/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                       MISCELLANEOUS STUFF                            %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/** converts an oligo into a unique integer; this function only works up to a
 * 14mer after that it will create an integer larger than will fit in 32bits
 * the function accepts a start so that a character buffer can be passed
 * when processing large sequences
 * only works for A,C,G,T  returns maximum value for any sequence containing
 * other IUPAC nucliec acid codes
 */
inline int BIO_seqToArrayIndex(const char *sequence, int start, int length);

/** pass the index for the previous word (from seqToArrayIndex), as well as
 * the next character, to get code for the new word. Used when moving 
 * along a sequence to avoid calculating the index for the whole word everytime
 */
inline int BIO_addCodeOne(int code, char c);

/** convert unique index from seqToArrayIndex into sequence */
inline void BIO_arrayIndexToSeq(int index, char *sequence, int MER_SIZE);

/* not for outside use */
int* getAAToINT();


/** counts the characters in any string (DNA, AA, text, etc...)
 * results are in global variable BIO_charCounts[]
 * to get the count of any character simply use BIO_charCounts['A'], 
 * BIO_charCounts['C'], etc...
 * if append is true then the previous counts are not cleared (so you can,
 * for example, obtain global counts for an entire alignment
 */
void BIO_getCounts(const char *s, int length, int append);

/** returns the translation of \a seq using the translation table 
 * corresponding to \a transCode, \a frame determines the reading frame;
 * 0 = no offset, 1 = offset by one, 2 = offset by two.
 * Note: returned char array MUST BE FREED.  To get all six translations you
 * need to call BIO_reverseComplement() and translate the reverse complement in 3 frames as well.
 */
char* BIO_getTranslation(const char *seq, BIO_geneticCode genCode, int frame);

/** returns an array which all ascii characters can be indexed into.  The values of the index are A = 0, C = 1, G = 2, T/U = 3, and the values are the same for any case (i.e. a==A).  All other characters are = -1.
*/
int* BIO_getNucToInt(void);


/** returns the melting temperature for the sequence; currently uses a primitive algorithm,
    but pairing energies should be added later.
   */
inline double BIO_getMeltingTemp(const char *seq, int start, int length);



/** checks various things to make sure the sequences are aligned 
 * (e.g. they must all be the same length) 
 * returns TRUE if it is an alignment
 */
inline int BIO_isAlignment(BIO_sequences s);

/** returns 1 if \a s contains only ACTG (case-insensitive) */
int BIO_isDNA(BIO_seqItem s);


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                FUNCTIONS FOR printing sequences (most moved to BIO_parsers.h  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/** prints S to a FASTA_NUM file num lines are of length \a lineWidth. 
 *  You must have an array of floating point numbers in your seqItem data structures
 *  for this to work. */
void BIO_printFASTA_NUM(FILE *outfile, BIO_sequences S, int lineWidth);

/** generates random nucleotide from uniform distribution */
inline char BIO_randNuc();

/** returns the length of the longest sequence in \a s */
inline int BIO_maxSeqLength(const BIO_sequences s);

/** gets the nucleotide at position \a row \a column in a matrix of sequences \a s 
 * if r or c is greater than the bounds of the sequence it returns a blank ' '
 */
inline char BIO_getNuc(BIO_sequences s, int row, int column);

/** returns the number of sequences in \a S */
inline int BIO_getNumSeqs(const BIO_sequences S);

/** gets the character array at \a row 
 * if r is greater than BIO_getNumSeqs(\a S) - 1 or less than 0 it returns NULL */
inline char* BIO_getSeq(BIO_sequences s, int row);

/** returns whether the sequence is nucleotide or amino acids */
inline BIO_seqType BIO_getSeqType(BIO_sequences s);

/** returns the name associated with the sequence at \a row in a  matrix
 * of sequences \a s; if \a s doesn't exist it returns NULL */
inline char* BIO_getSeqName(BIO_sequences s, int row);

/** returns the name of \a sI */
inline char* BIO_getSeqItemName(const BIO_seqItem sI);

/** sets the name of \a sI to \a newName; also updates the short name  */
void BIO_setSeqItemName(BIO_seqItem sI, const char *newName);

/** returns the short name of \a sI */
inline char* BIO_getSeqItemShortName(const BIO_seqItem sI);

/** returns the length of the sequence in \a sI */
inline int BIO_getSeqItemLength(const BIO_seqItem sI);

/** returns the sequence in \a sI */
inline char* BIO_getSeqItemSequence(const BIO_seqItem sI);

/** returns true if c is the complement of d (eg G-C is true G-A is false)
 * Rna complement regards G-U as true */
inline int BIO_isRnaComplement(char c, char d);

/** returns true if c is the complement of d (eg G-C is true G-A is false) */
inline int BIO_isDnaComplement(char c, char d);


/** returns the abbreviated name associated with the sequence at \a row \a column in a matrix
 * of sequences \a s; if sequence doesn't exist it returns NULL */
inline char* BIO_getShortSeqName(BIO_sequences S, int row);

/** sets the length to use for the shortened sequence name but ONLY in \a S. 
    It also updates all previous names.  To change the default short sequence name for all
    future sequences use BIO_setDefaultShortNameLength() */
void BIO_setShortNameLength(BIO_sequences S, int length);

/** sets the length to use for the shortened sequence name.  This number can also be set
    locally for a specific BIO_sequences structure using BIO_setShortNameLength() */
void BIO_setDefaultShortNameLength(int length); 

/** sets the user defined amino acid translation table so that calling 
 * translation routines with BIO_USER_DEFINED_CODE uses \a code.  
 * \a code should be an array of characters of the form: <pre>
char MY_TRANSLATION_TABLE[64] = {
        AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT,
        AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT,

        CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, 
        CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT,

        GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT,
        GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT,

	TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT,
	TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
};
</pre>
but replacing the codons with their corresponding amino acids.

Note: be careful when allocating \a code if you create the array
in a function and the memory removed when out of scope, it will create
a difficult to find bug/ program crash if you then translate something with
BIO_USER_DEFINED_CODE
 */
void BIO_setUserDefinedCode(char code[]);

/** works only for an alignment matrix;
 * returned array has a 1 where the column contains more than \a frequency gaps
 *(BIO_ANY_GAP to mark every column that has a gap; use BIO_ALL_GAPPED
 * to mark only columns that are all gaps)
 */
int* BIO_getGapSites(BIO_sequences s, double frequency);


/** returns a \a type consensus sequence from the alignment in \a s. If an error occurs
 * (e.g. if the sequences aren't aligned), it returns \a NULL. The function excludes (puts a gap for) 
 * any column where \a exclude[i] is true; this is useful for defining gapped columns
 * using BIO_getGapSites() or your own function. Pass \a NULL to attempt creation of a 
 * consensus character for all sites.
 */
char* BIO_getConsensus(BIO_sequences s, BIO_consensusType type, const int *exclude);


__inline__ char* createShortName(const char *longName, int nameLength);


#ifdef __cplusplus
}
#endif

#endif // BIO_SEQUENCE_H
