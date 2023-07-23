#include "BIO_sequence.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#define MAX(a,b) (a > b) ? a : b

/* macro that calculates melting temperature using the Wallace method 
 * Wallace RB, Shaffer J, Murphy RF, Bonner J, Hirose T, and Itakura, K. 
 * Nucleic Acids Research 6, 3543 (1979) */
#define WALLACE(A,B,C,D) 2.0*((A)+(B)) + 4.0*((C)+(D))

int DEFAULT_SHORT_NAME_LENGTH = 10;

double wallaceTemp(const char *seq, int offset, int length);
void getMatrix(BIO_sequences S);
/* infers whether sequence is amino acids or nucleotide using some rudimentary statistics */
void inferSeqType(BIO_sequences s);
/* a global reference to the current type for a gap; use this so in the future the program
 * is not limited to using '-' as a gap character */
int BIO_GAP = '-';



BIO_seqItem AppendSequence(BIO_seqItem headRef, BIO_seqItem newNode);
void printList(BIO_seqItem mySeqs);

__inline__ void clearSeqCounts();
int BIO_seqBufferSize=BIO_SEQ_BUFFER_SIZE;

/* returns the consensus character for <sType> according to <conType> using
 * <charCounts> */
__inline__ char getConsensusChar(BIO_seqType sType, BIO_consensusType conType, int *charCounts); 
/* returns the nucleotide character for <conType> using <charCounts> */
__inline__ char getNucConsensusChar(BIO_consensusType conType, const int *charCounts);
/* returns the amino acid character for <conType> using <charCounts> */
__inline__ char getAAConsensusChar(BIO_consensusType conType, const int *cc);

/* updates all the short names for all the sequences */
void updateShortNames(BIO_sequences S, int length);

/* initialize global variable used when speed is needed */
int merAND=0;
char NUCLEOTIDES[4]={'A','C','G','T'};

int BIO_NUM_NEGATIVE_CHARGE = 2;
int BIO_NUM_POSITIVE_CHARGE = 3;
int BIO_NUM_CHARGED         = 5;
int BIO_NUM_NONPOLAR        = 7;
int BIO_NUM_POLAR           = 8;

int BIO_charCounts[SEQ_NUM_CHAR_COUNTS];

#define ROW_GAP 100
#define COL_GAP 101
#define NO_GAP  102
#define END_ALN 103


/* best for detecting long and weak alignments */
/* from the durbin, eddy, krogh, mitchison book */
int BLOSUM_50[20][20] = 
{       /*  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V */ 
	{   5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0  }, /* A */
	{  -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3  }, /* R */
	{  -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3  }, /* N */
	{  -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4  }, /* D */
	{  -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1  }, /* C */
	{  -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3  }, /* Q */
	{  -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3  }, /* E */
	{   0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4  }, /* G */
	{  -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4  }, /* H */
	{  -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1, -3  }, /* I */
	{  -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1  }, /* L */
	{  -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3  }, /* K */
	{  -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1  }, /* M */
	{  -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1  }, /* F */
	{  -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3  }, /* P */
	{   1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2  }, /* S */
	{   0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0  }, /* T */
	{  -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3  }, /* W */
	{  -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1  }, /* Y */
	{   0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5  }, /* V */
};





/* this array allows you to use a nucleotide as an array value to obtain
 * its integer equivalent:  A=a==0, C=c=1, G=g=2, T=t=U=u=3
 */
static int NUC_TO_INT[255]={
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

/* amino acids corresponding to codon indexes obtained by seqToArrayIndex for BIO_STANDARD_CODE */
char STANDARD_INT_TO_AA[64] = {
	'K', /* AAA */ 'N', /* AAC */ 'K', /* AAG */ 'N', /* AAT */
	'T', /* ACA */ 'T', /* ACC */ 'T', /* ACG */ 'T', /* ACT */
	'R', /* AGA */ 'S', /* AGC */ 'R', /* AGG */ 'S', /* AGT */
	'I', /* ATA */ 'I', /* ATC */ 'M', /* ATG */ 'I', /* ATT */

	'Q', /* CAA */ 'H', /* CAC */ 'Q', /* CAG */ 'H', /* CAT */
	'P', /* CCA */ 'P', /* CCC */ 'P', /* CCG */ 'P', /* CCT */
	'R', /* CGA */ 'R', /* CGC */ 'R', /* CGG */ 'R', /* CGT */
	'L', /* CTA */ 'L', /* CTC */ 'L', /* CTG */ 'L', /* CTT */

	'E', /* GAA */ 'D', /* GAC */ 'E', /* GAG */ 'D', /* GAT */
	'A', /* GCA */ 'A', /* GCC */ 'A', /* GCG */ 'A', /* GCT */
	'G', /* GGA */ 'G', /* GGC */ 'G', /* GGG */ 'G', /* GGT */
	'V', /* GTA */ 'V', /* GTC */ 'V', /* GTG */ 'V', /* GTT */

	'*', /* TAA */ 'Y', /* TAC */ '*', /* TAG */ 'Y', /* TAT */
	'S', /* TCA */ 'S', /* TCC */ 'S', /* TCG */ 'S', /* TCT */
	'*', /* TGA */ 'C', /* TGC */ 'W', /* TGG */ 'C', /* TGT */
	'L', /* TTA */ 'F', /* TTC */ 'L', /* TTG */ 'F', /* TTT */
};

/* amino acids corresponding to codon indexes obtained by seqToArrayIndex for BIO_VERTEBRATE_MT_CODE */
char VERTEBRATE_MT_INT_TO_AA[64] = {
	'K', /* AAA */ 'N', /* AAC */ 'K', /* AAG */ 'N', /* AAT */
	'T', /* ACA */ 'T', /* ACC */ 'T', /* ACG */ 'T', /* ACT */
	'*', /* AGA */ 'S', /* AGC */ '*', /* AGG */ 'S', /* AGT */
	'M', /* ATA */ 'I', /* ATC */ 'M', /* ATG */ 'I', /* ATT */

	'Q', /* CAA */ 'H', /* CAC */ 'Q', /* CAG */ 'H', /* CAT */
	'P', /* CCA */ 'P', /* CCC */ 'P', /* CCG */ 'P', /* CCT */
	'R', /* CGA */ 'R', /* CGC */ 'R', /* CGG */ 'R', /* CGT */
	'L', /* CTA */ 'L', /* CTC */ 'L', /* CTG */ 'L', /* CTT */

	'E', /* GAA */ 'D', /* GAC */ 'E', /* GAG */ 'D', /* GAT */
	'A', /* GCA */ 'A', /* GCC */ 'A', /* GCG */ 'A', /* GCT */
	'G', /* GGA */ 'G', /* GGC */ 'G', /* GGG */ 'G', /* GGT */
	'V', /* GTA */ 'V', /* GTC */ 'V', /* GTG */ 'V', /* GTT */

	'*', /* TAA */ 'Y', /* TAC */ '*', /* TAG */ 'Y', /* TAT */
	'S', /* TCA */ 'S', /* TCC */ 'S', /* TCG */ 'S', /* TCT */
	'W', /* TGA */ 'C', /* TGC */ 'W', /* TGG */ 'C', /* TGT */
	'L', /* TTA */ 'F', /* TTC */ 'L', /* TTG */ 'F', /* TTT */
};

/* amino acids corresponding to codon indexes obtained by seqToArrayIndex for BIO_YEAST_MT_CODE 
 */
char YEAST_MT_INT_TO_AA[64] = {
	'K', /* AAA */ 'N', /* AAC */ 'K', /* AAG */ 'N', /* AAT */
	'T', /* ACA */ 'T', /* ACC */ 'T', /* ACG */ 'T', /* ACT */
	'R', /* AGA */ 'S', /* AGC */ 'R', /* AGG */ 'S', /* AGT */
	'M', /* ATA */ 'I', /* ATC */ 'M', /* ATG */ 'I', /* ATT */

	'Q', /* CAA */ 'H', /* CAC */ 'Q', /* CAG */ 'H', /* CAT */
	'P', /* CCA */ 'P', /* CCC */ 'P', /* CCG */ 'P', /* CCT */
	'X', /* CGA */ 'X', /* CGC */ 'R', /* CGG */ 'R', /* CGT */
	'T', /* CTA */ 'T', /* CTC */ 'T', /* CTG */ 'T', /* CTT */

	'E', /* GAA */ 'D', /* GAC */ 'E', /* GAG */ 'D', /* GAT */
	'A', /* GCA */ 'A', /* GCC */ 'A', /* GCG */ 'A', /* GCT */
	'G', /* GGA */ 'G', /* GGC */ 'G', /* GGG */ 'G', /* GGT */
	'V', /* GTA */ 'V', /* GTC */ 'V', /* GTG */ 'V', /* GTT */

	'*', /* TAA */ 'Y', /* TAC */ '*', /* TAG */ 'Y', /* TAT */
	'S', /* TCA */ 'S', /* TCC */ 'S', /* TCG */ 'S', /* TCT */
	'W', /* TGA */ 'C', /* TGC */ 'W', /* TGG */ 'C', /* TGT */
	'L', /* TTA */ 'F', /* TTC */ 'L', /* TTG */ 'F', /* TTT */
};
/* amino acids corresponding to codon indexes obtained by seqToArrayIndex for BIO_INVERTEBRATE_MT_CODE 
 */
char INVERTEBRATE_MT_INT_TO_AA[64] = {
	'K', /* AAA */ 'N', /* AAC */ 'K', /* AAG */ 'N', /* AAT */
	'T', /* ACA */ 'T', /* ACC */ 'T', /* ACG */ 'T', /* ACT */
	'S', /* AGA */ 'S', /* AGC */ 'S', /* AGG */ 'S', /* AGT */
	'M', /* ATA */ 'I', /* ATC */ 'M', /* ATG */ 'I', /* ATT */

	'Q', /* CAA */ 'H', /* CAC */ 'Q', /* CAG */ 'H', /* CAT */
	'P', /* CCA */ 'P', /* CCC */ 'P', /* CCG */ 'P', /* CCT */
	'R', /* CGA */ 'R', /* CGC */ 'R', /* CGG */ 'R', /* CGT */
	'L', /* CTA */ 'L', /* CTC */ 'L', /* CTG */ 'L', /* CTT */

	'E', /* GAA */ 'D', /* GAC */ 'E', /* GAG */ 'D', /* GAT */
	'A', /* GCA */ 'A', /* GCC */ 'A', /* GCG */ 'A', /* GCT */
	'G', /* GGA */ 'G', /* GGC */ 'G', /* GGG */ 'G', /* GGT */
	'V', /* GTA */ 'V', /* GTC */ 'V', /* GTG */ 'V', /* GTT */

	'*', /* TAA */ 'Y', /* TAC */ '*', /* TAG */ 'Y', /* TAT */
	'S', /* TCA */ 'S', /* TCC */ 'S', /* TCG */ 'S', /* TCT */
	'W', /* TGA */ 'C', /* TGC */ 'W', /* TGG */ 'C', /* TGT */
	'L', /* TTA */ 'F', /* TTC */ 'L', /* TTG */ 'F', /* TTT */
};

/* allow user to create their own translation table */
char *USER_DEFINED_INT_TO_AA = NULL;

/* IUPAC nucleotide complements */
char COMPLEMENT[255]={
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,'-','.',-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,'T','V','G','H',-1,-1,'C','D',-1,-1,'.',-1,'K','N',-1,-1,-1,'Y','S','A','A','B','W','X','R',
	-1,-1,-1,-1,'^',-1,-1,'t','v','g','h',-1,-1,'c','d',-1,-1,'m',-1,'k','n',-1,-1,-1,'y','s','a','a','b','w',
	'x','r',-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};


/* reverses string s in place;  this is super cool, unfortunately its
 * not original; algorithm comes form the K&R book
 * note this will crash program if someone has used char *myString="ACTACTACGA"
 * to initialize the array; must use malloc or char myString[]="ACTACTACGA"
 */
void BIO_reverseString(char s[]) {
	int i,j,c;

	for(i=0, j=strlen(s)-1; i<j; i++, j--)
		c=s[i], s[i]=s[j], s[j]=c;
}

void BIO_stringToUpper(char *string) {
	int i;
	int length = strlen(string);

	for (i=0; i<length; i++)
		string[i] = toupper(string[i]);
}

void BIO_stringToLower(char *string) {
	int i;
	int length = strlen(string);

	for (i=0; i<length; i++)
		string[i] = tolower(string[i]);
}

void BIO_reverseComplement(char s[]) {
	int length,i;
	length=strlen(s);
	//printf("length is %d\n", length);
	for (i=0; i<length; i++)
		s[i]=COMPLEMENT[(int)s[i]];

	// before this was commented out for some reason? 
	BIO_reverseString(s);
}

char* BIO_getTranslation(const char *seq, BIO_geneticCode transCode, int frame) {
	int i;
	int length = strlen(seq);
	char *trans;
	char *INT_TO_AA=NULL;
	char codon[3];
	int score;
	int nucCount=0;

	switch(transCode) {
		case BIO_STANDARD_CODE:          INT_TO_AA = STANDARD_INT_TO_AA;
			break;
		case BIO_VERTEBRATE_MT_CODE:     INT_TO_AA = VERTEBRATE_MT_INT_TO_AA;
			break;
		case BIO_INVERTEBRATE_MT_CODE:   INT_TO_AA = INVERTEBRATE_MT_INT_TO_AA;
			break;
		case BIO_YEAST_MT_CODE:          INT_TO_AA = YEAST_MT_INT_TO_AA;
			break;
		case BIO_USER_DEFINED_CODE:      
			if (USER_DEFINED_INT_TO_AA != NULL)
				INT_TO_AA = USER_DEFINED_INT_TO_AA;
			else /* use the standard code if the user didn't give us one */
				INT_TO_AA = STANDARD_INT_TO_AA;
			break;
		default:                         INT_TO_AA = STANDARD_INT_TO_AA;
			break;
	}

	if (frame < 0 || frame > 2) {
		fprintf(stderr, "frame must be between 0 and 2\n");
		return NULL;
	}

	trans=(char*)malloc( sizeof(char) * ( ( (length-frame)/3) + 1 ));
	for (i=frame; i<length; i+=3) {
		codon[0]=seq[i];
		codon[1]=seq[i+1];
		codon[2]=seq[i+2];
		score = BIO_seqToArrayIndex(codon, 0, 3);
		if (score < 64 && score >= 0) {
		//	printf("score is %2d codon is %c%c%c AA is %c\n", score, codon[0],codon[1],codon[2],INT_TO_AA[score]);	
			trans[nucCount++] = INT_TO_AA[score];
		}
	}
	trans[nucCount] = '\0';
	return trans;
}




/* this is slow, but it works */
BIO_seqItem BIO_iupacToACGT(BIO_seqItem mySeqs) {
	int i;
	BIO_seqItem one,two,three,four,temp;
	char *seq;
	int length;
	int finished=0;
	int last=0; // used to prevent looking at nucleotide we know are fixed

	// compute length outside loop
	length=strlen(mySeqs->seq);

	while (!finished) {
		temp=mySeqs;
		mySeqs=temp->next;
	
		seq=temp->seq;
		
		for (i=last; i<=length; i++) {
			last=i;	// remember where we left off to start there next time

			/* not every conditional 'breaks' because we only fix one
			 * nucleotide at a time */

			/* split puRine into Guanine and Adenine */
			if (seq[i] == 'R') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='A';
				two->seq[i]='G';

				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);

				free(seq);
				free(temp);
				break;
			}
			/* split pYrimidine into Guanine and Adenine */
			else if (seq[i] == 'Y') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='C';
				two->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				free(seq);
				free(temp);
				break;
			}
			/* split aMino into Adenine and Cytosine */
			else if (seq[i] == 'M') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='A';
				two->seq[i]='C';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);

				free(seq);
				free(temp);
				break;
			}
			/* split Keto into Guanine and Tyrosine */
			else if (seq[i] == 'K') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='G';
				two->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);

				free(seq);
				free(temp);
				break;
			}
			/* split Strong (3 H bonds) into Guanine and Cytosine */
			else if (seq[i] == 'S') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='G';
				two->seq[i]='C';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);

				free(seq);
				free(temp);
				break;
			}
			/* split Weak (2 H bonds) into Adenine and Tyrosine */
			else if (seq[i] == 'W') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				one->seq=strdup(seq);
				two->seq=strdup(seq);
				one->seq[i]='A';
				two->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);

				free(seq);
				free(temp);
				break;
			}
			/* split H (not-G) A,C,T */
			else if (seq[i] == 'H') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				three=BIO_initSeqItem();
				one->seq=strdup(seq), two->seq=strdup(seq), three->seq=strdup(seq);
				one->seq[i]='A', two->seq[i]='C', three->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				mySeqs=AppendSequence(mySeqs,three);

				free(seq);
				free(temp);
				break;
			}
			/* split B (not-A) C,G,T */
			else if (seq[i] == 'B') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				three=BIO_initSeqItem();
				one->seq=strdup(seq), two->seq=strdup(seq), three->seq=strdup(seq);
				one->seq[i]='C', two->seq[i]='G', three->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				mySeqs=AppendSequence(mySeqs,three);

				free(seq);
				free(temp);
				break;
			}
			/* split V (not-U) G,C,A */
			else if (seq[i] == 'V') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				three=BIO_initSeqItem();
				one->seq=strdup(seq);   two->seq=strdup(seq);   three->seq=strdup(seq);
				one->seq[i]='A'; two->seq[i]='C'; three->seq[i]='G';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				mySeqs=AppendSequence(mySeqs,three);

				free(seq);
				free(temp);
				break;
			}
			/* split D (not-C) G,A,T */
			else if (seq[i] == 'D') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				three=BIO_initSeqItem();
				one->seq=strdup(seq);   two->seq=strdup(seq);   three->seq=strdup(seq);
				one->seq[i]='A'; two->seq[i]='G'; three->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				mySeqs=AppendSequence(mySeqs,three);

				free(seq);
				free(temp);
				break;
			}
			/* split N (aNy) A,C,G,T */
			else if (seq[i] == 'N') {
				// allocate
				one=BIO_initSeqItem();
				two=BIO_initSeqItem();
				three=BIO_initSeqItem();
				four =BIO_initSeqItem();
				one->seq=strdup(seq);   two->seq=strdup(seq);   
				three->seq=strdup(seq); four->seq=strdup(seq);
				one->seq[i]='A';   two->seq[i]='C'; 
				three->seq[i]='G'; four->seq[i]='T';
				mySeqs=AppendSequence(mySeqs,one);
				mySeqs=AppendSequence(mySeqs,two);
				mySeqs=AppendSequence(mySeqs,three);
				mySeqs=AppendSequence(mySeqs,four);

				free(seq);
				free(temp);
				break;
			}
			// if we make it to the last character we are done
			else if (seq[i]=='\0')
				finished=1;

		}
	}

	return temp;
}

BIO_seqItem BIO_initSeqItem() {
	BIO_seqItem sI;
	sI=(BIO_seqItem)malloc(sizeof(*sI));
	sI->seq=NULL;
	sI->name=NULL;
	sI->shortName=NULL;
	sI->next=NULL;
	sI->DATA=NULL;
	sI->seqQual=NULL;
	sI->bufferSize=0;
	sI->seqLength=0;

	return sI;
}

BIO_sequences BIO_initSequences() {
	BIO_sequences s = (BIO_sequences)malloc(sizeof(*s));
	s->seqs=NULL;	
	s->numSeqs=0;
	s->hasMatrix=0;
	s->shortNameLength=10;
	s->matrix=NULL;
	s->seqType=BIO_UNCHECKED_TYPE;

	return s;
}

BIO_sequences BIO_copySequences(const BIO_sequences old) {
	BIO_seqItem current;

	BIO_sequences s = (BIO_sequences)malloc(sizeof(*s));
	s->seqs            = old->seqs;
	s->numSeqs         = old->numSeqs;
	s->hasMatrix       = old->hasMatrix;
	s->shortNameLength = old->shortNameLength;
	s->matrix          = NULL;
	s->seqType         = old->seqType;

	for (current=s->seqs; current != NULL; current=current->next)
		BIO_addSeqItem(s, BIO_copySeqItem(current));

	return s;
}

BIO_seqItem BIO_copySeqItem(const BIO_seqItem old) {
	BIO_seqItem sI = NULL;
	if (!old) return sI;

	sI=(BIO_seqItem)malloc(sizeof(*sI));
	sI->seq        = strdup(old->seq);
	sI->name       = strdup(old->name);
	sI->shortName  = strdup(old->shortName);
	sI->next       = NULL; /* dangerous to mess with this */
	sI->DATA       = NULL; /* don't know how to copy this */
	sI->bufferSize = old->bufferSize;
	sI->seqLength  = old->seqLength;

	return sI;
}

void BIO_setShortNameLength(BIO_sequences S, int length) { 
	if (S->shortNameLength !=  length)
		updateShortNames(S, length);

	S->shortNameLength=length; 
}

void BIO_setUserDefinedCode(char code[]) {
	USER_DEFINED_INT_TO_AA = code;
}

BIO_seqItem BIO_getSeqItem(BIO_sequences S, const char *name, int BOOLEAN) {
	BIO_seqItem current;
	for (current=S->seqs; current != NULL; current=current->next)
		if(strcmp(current->name, name) == 0)
			return current;
 
	if (BOOLEAN) {
		S->numSeqs++;
		current = BIO_addSeqItem(S, BIO_createSeqItem(name, NULL));
		/* give it a name */

		if (BIO_seqBufferSize) {
			current->seq = (char *)malloc( sizeof(char) * BIO_seqBufferSize);
			current->bufferSize = BIO_seqBufferSize;
		}
		return current;
	}
	return NULL;
}

BIO_seqItem BIO_createSeqItem(const char *name, const char *seq) {
	BIO_seqItem newS = NULL;

	if (!name)
		return NULL;

	newS = BIO_initSeqItem();
	if (seq!=NULL) {
		newS->seq       = strdup(seq);
		newS->seqLength = strlen(newS->seq);;
		newS->bufferSize = newS->seqLength+1;
	}

	newS->name = strdup(name);
	newS->shortName = createShortName(name, DEFAULT_SHORT_NAME_LENGTH);

	return newS;
}

void BIO_setDefaultShortNameLength(int length) 
{ DEFAULT_SHORT_NAME_LENGTH = length; }

void BIO_catStr(BIO_seqItem sL, const char *seq, int length) {
	int i=0;
	int oldEnd = sL->seqLength;

	if ((length + sL->seqLength) >= sL->bufferSize) {
		sL->bufferSize = (sL->bufferSize + length) * 2;
		sL->seq = (char*)realloc(sL->seq, sizeof(char) * sL->bufferSize);
//		printf("buffer size is %d\n", sL->bufferSize);
	}

	/* concatenate if there was sequence before */
	if (sL->seqLength) {
		/* strcat get too slow with giant sequences like chromosomes 
		   because I think strcat has to scan the entire sequence for the '\0'
	           but since we know where the end is, we can simplify the process */
		//	strcat(sL->seq, seq);
//		printf("first %s\n", sL->seq);
//		printf("second %s\n", seq);
		i=0;
		while (seq[i] != '\0')
			sL->seq[oldEnd++] = seq[i++];

		sL->seq[oldEnd] = '\0';
//		printf("after %s\n", sL->seq);
	}
	/* copy if it is a new sequence */
	else {
		strcpy(sL->seq, seq);
	}
	
	sL->seqLength += length;
}

void BIO_seqClean(BIO_sequences S) {
	BIO_seqItem current;
	char *temp;
	for (current=S->seqs; current!=NULL; current=current->next) {
		if (current->seqLength+1 < current->bufferSize) {
			temp=current->seq;
			current->seq = (char*)malloc(sizeof(char) * (current->seqLength+1) );
			strcpy(current->seq,temp);
			free(temp);
		}
	}
}


void BIO_destroySequences(BIO_sequences S) {
	BIO_seqItem current;
	BIO_seqItem temp;
	for (current=S->seqs; current!=NULL; current=temp) {
		temp=current->next;
		BIO_destroySeqItem(current);
	}

	if (S->hasMatrix)
		free(S->matrix);

	free(S);
}

void BIO_destroySeqItem(BIO_seqItem sI) {
	if (sI == NULL) 
		return;
	free(sI->seq);
	free(sI->seqQual);
	free(sI->name);
	free(sI->shortName);
	free(sI->DATA);
	free(sI);
}



BIO_seqItem BIO_addSeqItem(BIO_sequences s, BIO_seqItem sI) {
	s->seqs = AppendSequence(s->seqs, sI);	

	if (strlen(sI->shortName) != s->shortNameLength) {
		free(sI->shortName);
		sI->shortName = createShortName(sI->name, s->shortNameLength);
	}

	return sI;
}

char* createShortName(const char *longName, int length)
{
	if (!length) {
		length = DEFAULT_SHORT_NAME_LENGTH;
	}
	char *shortName = (char*)malloc(sizeof(char) * (length + 1));	
	strncpy(shortName, longName, length);	
	shortName[length]='\0';
	return shortName;
}

void updateShortNames(BIO_sequences S, int length) {
	BIO_seqItem current;

	for (current=S->seqs; current!=NULL; current=current->next) {
		free(current->shortName);
		current->shortName = createShortName(current->name, length);
	}
}

BIO_seqItem AppendSequence(BIO_seqItem headRef, BIO_seqItem newNode) {
	BIO_seqItem current;
	current=headRef;	

	if (current == NULL) 
		headRef=newNode;
	else {
		while(current->next != NULL)
			current=current->next;

		current->next = newNode;
	}

	return headRef;
}


void BIO_printFASTA_NUM(FILE *out, BIO_sequences S, int width) {
	BIO_seqItem current;
	float *currentNums;
	int i;
	int count=0;
	for (current=S->seqs; current!=NULL; current=current->next) {
		count=0;
		if (current->DATA != NULL) {
			fprintf(out,">%s\n", current->name); 
			currentNums = (float*)current->DATA;
			for (i=0; i<current->seqLength; i++) {
				fprintf(out, "%f", currentNums[i]);
				if (++count % width == 0 && (i!=current->seqLength-1))
					putc('\n', out);
				else 
					putc(' ', out);
			}
			putc('\n', out);
			putc('\n', out);
		}
	}
}



void printList(BIO_seqItem mySeqs) {
		BIO_seqItem current;
		printf("LIST\n");
		for (current=mySeqs; current!=NULL; current=current->next)
			printf("S => sequence is %s\n", current->seq);
}

//void arrayIndexToSeq

/* convert nuc sequence into a unique index */
int BIO_seqToArrayIndex(const char *seq, int start, int length) {
	int i;
	int stop=start+length;
	int index=0;
//	extern int NUC_TO_INT[];

	for (i=start; i<stop; i++) {
		//putchar(seq[i]);
		index *= 4; /* shifts index two bits (just the right size for two nucs */
		// stop if we see nucleotide that is not ACGT
		if (NUC_TO_INT[(int)seq[i]] == -1)
			return -1;
	
		index += NUC_TO_INT[(int)seq[i]];
	}
	return index;
}


int BIO_addCodeOne(int code, char c) {
	//extern int NUC_TO_INT[];
	extern int merAND;
	int C=(int)c;

	// return -1 for rubbish nucleotides
	if (NUC_TO_INT[C] == -1)
		return -1;

	code <<= 2;	// like multiplication by 4 but faster
	code=(code & merAND);
	code += NUC_TO_INT[C];
	return code;
}

/* works for AA and nucs just counts characters in string */
void BIO_getCounts(const char *s, int length, int append) {
	int i;
	//extern int BIO_charCounts[];
	if (!append)
		clearSeqCounts();

	for (i=0; i<length; i++) {
		BIO_charCounts[ (int)s[i] ]++;
	}
}

void clearSeqCounts() {
	int i;
	//extern int BIO_charCounts[];
	for (i=0; i<SEQ_NUM_CHAR_COUNTS; i++)
		BIO_charCounts[i]=0;
}

/* convert unique index from seqToArrayIndex into sequence */
void BIO_arrayIndexToSeq(int index, char *seq, int MER) {
	int i,j;
	//extern char NUCLEOTIDES[];
	for (j=0; j<MER; j++) {
		i=(index & 3);
		seq[MER-j-1]=NUCLEOTIDES[i];
		index=index>>2;
	}
	seq[MER]='\0';
}

int factorial(int n)
{
	int i=0;
	int product=1;

	if (n<0)
		fprintf(stderr, "\nError: factorial of negative integer not defined\n");

	for (i=n; i>1; --i)
		product=product*i;

	return product;
}

double BIO_getMeltingTemp(const char *seq, int start, int length) {
	double temp=0.0;
	//if (length <= 20)
		return wallaceTemp(seq, start, length);

	return temp;
}

int* BIO_getNucToInt(void) {
	return NUC_TO_INT;
}

double wallaceTemp(const char *seq, int offset, int length) {
	int i;
	char c;
	double A,C,G,T;
	A=C=G=T=0.0;
	for (i=0; i<length; i++) {
		c=seq[i+offset];
		/* quit if we hit the end of the string */
		if (c == '\0')
			//return WALLACE(A,C,G,T);
			return WALLACE(A,C,G,T);

		switch(c) {
			case 'A': A++; break;
			case 'C': C++; break;
			case 'G': G++; break;
			case 'T': G++; break;
			case 'R': G+= 0.5, A += 0.5; break;
			case 'Y': T+= 0.5, C += 0.5; break;
			case 'M': A+= 0.5, C += 0.5; break;
			case 'K': G+= 0.5, T += 0.5; break;
			case 'S': G+= 0.5, C += 0.5; break;
			case 'W': A+= 0.5, T += 0.5; break;
			case 'H': A+= 1.0/3.0, C += 1.0/3.0, T += 1.0/3.0; break;
			case 'B': G+= 1.0/3.0, C += 1.0/3.0, T += 1.0/3.0; break;
			case 'V': G+= 1.0/3.0, C += 1.0/3.0, A += 1.0/3.0; break;
			case 'D': G+= 1.0/3.0, A += 1.0/3.0, T += 1.0/3.0; break;
			case 'N': A += 0.25, C += 0.25, G += 0.25, T+=0.25; break;
		        /* not sure what to do if I don't understand the letter
			 * is it a gap? amino acids? just junk?  leaving it blank for now */
			default: break;
                }
        }
	return WALLACE(A,C,G,T);
}

char BIO_randNuc() 
{ return NUCLEOTIDES[ rand() % 4 ] ; } 

char BIO_getNuc(BIO_sequences s, int r, int c) {
	if (!s->hasMatrix)
		getMatrix(s);
	
	if (r >=0 && r < s->numSeqs && c >=0 && c < s->matrix[r]->seqLength)
		return s->matrix[r]->seq[c];
	else 
		return ' ';
}

char* BIO_getSeq(BIO_sequences s, int r) {
	if (!s->hasMatrix)
		getMatrix(s);

	if (r >=0 && r < s->numSeqs)
		return s->matrix[r]->seq;

	return NULL;
}

BIO_seqType BIO_getSeqType(BIO_sequences s) {
	if (!s->seqType)
		inferSeqType(s);

	return s->seqType;
}

void inferSeqType(BIO_sequences S) {
	BIO_seqItem current;
	int numNucs=0;
	int ACGT=0;

	clearSeqCounts();

	for (current=S->seqs; current!=NULL && numNucs < 10000; current=current->next) {
		BIO_getCounts(current->seq, current->seqLength, 1);
		numNucs += current->seqLength;
	}

	/* disregard gaps */
	numNucs -= BIO_charCounts['-'];
	numNucs -= BIO_charCounts['.'];
	numNucs -= BIO_charCounts['_'];

	ACGT += BIO_charCounts['A'];
	ACGT += BIO_charCounts['C'];
	ACGT += BIO_charCounts['G'];
	ACGT += BIO_charCounts['T'];


	/* if more than half of the characters are ACGT assume the sequence is nucleotides */
	if ( (ACGT/(double)numNucs) < 0.5 )
		S->seqType = BIO_AMINO_ACID;
	else
		S->seqType = BIO_NUCLEOTIDE;

}

BIO_seqItem BIO_getSeqItemByIndex(BIO_sequences s, int row) {
	if (!s->hasMatrix)
		getMatrix(s);

	if (row >=0 && row < s->numSeqs)
		return s->matrix[row];
	else
		return NULL;
}

char* BIO_getSeqName(BIO_sequences s, int r) {
	if (!s->hasMatrix)
		getMatrix(s);

	if (r >=0 && r < s->numSeqs)
		return BIO_getSeqItemName(s->matrix[r]);
	else
		return NULL;
}

char* BIO_getShortSeqName(BIO_sequences s, int r) {
	if (!s->hasMatrix)
		getMatrix(s);

	if (r >=0 && r < s->numSeqs) {
		//fprintf(stderr, "name %s nuc is %c\n", s->matrix[r]->name, s->matrix[r]->seq[c]);
		return BIO_getSeqItemShortName(s->matrix[r]);
	}
	else {
		return NULL;
	}
}

void getMatrix(BIO_sequences S) {
	BIO_seqItem current;
	int i=0;
	fprintf(stderr, "building matrix\n");

	if (S->matrix != NULL)
		free(S->matrix);

	S->matrix = (BIO_seqItem*)malloc(sizeof(*S->matrix) * S->numSeqs);	
	for (current=S->seqs; current!=NULL; current=current->next)
		S->matrix[i++]=current;

	S->hasMatrix=1;
}

int BIO_maxSeqLength(const BIO_sequences S) {
	BIO_seqItem current;
	int maxLength=0;
	fprintf(stderr, "starting bio %d\n", maxLength);
	for (current=S->seqs; current!=NULL; current=current->next) {
		if (current->seqLength > maxLength)
			maxLength=current->seqLength;
	}
	fprintf(stderr, "ending bio %d\n", maxLength);

	return maxLength;
}

int BIO_isAlignment(BIO_sequences S) {
	int length = S->seqs->seqLength;
	BIO_seqItem current;
	for (current=S->seqs; current!=NULL; current=current->next) 
		if (current->seqLength != length)
			return -1;

	return 1;
}

int* BIO_getGapSites(BIO_sequences S, double frequency) {
	int *gappedSites=NULL;
	int noMoreThan=0;
	double x=-999999.0;
	int numGaps,i,j;

	/* can't look for gapped columns if S isn't an alignment */
	if (!BIO_isAlignment(S))
		return NULL;
	/* can't look for a negative number of gaps */
	if (frequency < 0.0)
		return NULL;

	if (!S->hasMatrix) 
		getMatrix(S);


	//gappedSites = (int*)calloc( S->seqs->seqLength, sizeof(*gappedSites) );
	gappedSites = (int*)calloc( S->seqs->seqLength, sizeof(int) );

	if (frequency == BIO_ANY_GAP)
		noMoreThan = 0;
	else if (frequency == BIO_ALL_GAPPED)
		noMoreThan = S->numSeqs-1;
	else {
		while (x < frequency) 
			x = ++noMoreThan / (double)S->numSeqs;

		noMoreThan--;
	}

	

	for (i=0; i<S->seqs->seqLength; i++) {
		numGaps=0; // no gaps in this column yet
		for (j=0; j<S->numSeqs; j++) {
			if (S->matrix[j]->seq[i] == BIO_GAP) {
				if (++numGaps > noMoreThan) {
					gappedSites[i]=1; 
					break; // call column gapped when it is greater <percent>
				}
			}
		}
	}

	return gappedSites;
}

char* BIO_getConsensus(BIO_sequences s, BIO_consensusType type, const int *exclude) {
	int i,j;
	char *consensus=NULL;

	/* can't build a consensus sequence if the sequences aren't aligned */
	if (!BIO_isAlignment(s))
		return consensus;
	/* if the sequence type hasn't been determined, sort that out now */
	if (s->seqType == BIO_UNCHECKED_TYPE)
		inferSeqType(s);
	/* can't build a consensus if we don't know if they are amino acids or proteins */
	if (s->seqType == BIO_UNKNOWN_TYPE)
		return consensus;

	consensus = (char*) malloc(sizeof(char) * s->seqs->seqLength + 1);

	if (!s->hasMatrix)
		getMatrix(s);


	/* set the last character */
	consensus[s->seqs->seqLength]='\0';
	for (i=0; i < s->seqs->seqLength; i++) {
		/* excluded rows are just gaps */
		if (exclude != NULL && exclude[i]) {
			consensus[i] = '-';
		}
		/* getCounts for the whole column */
		else {
			clearSeqCounts();
			for (j=0; j < s->numSeqs; j++) 
				BIO_charCounts[(int)s->matrix[j]->seq[i]]++;

			consensus[i] = getConsensusChar(s->seqType, type, BIO_charCounts);
		}
	}

	return consensus;
}

char getConsensusChar(BIO_seqType sType, BIO_consensusType conType, int *charCounts) {
	switch(sType) {
		case BIO_NUCLEOTIDE:
			return getNucConsensusChar(conType, charCounts); break;

		case BIO_AMINO_ACID:
			return getAAConsensusChar(conType, charCounts); break;
	}

	return '-'; /* if nothing else return a gap */
}

char getAAConsensusChar(BIO_consensusType conType, const int *cc) {
	int i;
	char AA[] = "ARNDCQEGHILKMFPSTWYV";
	char U[] = "DE";
	char _O[] = "KHR";
	char _Z[] = "DEKRH";
	char B[] = "AIVLFPM";
	char _J[] = "STYHCNQW";
	char con = '-';
	int count=0;
	int localCount = 0;
	int sum = 0;
	/* I purposely set this to zero rather than negative
	 * to catch weird things in BIO_MAJORITY_RULE if no ACGT is
	 * present */
	int max = 0;

	switch(conType) {
		case BIO_STRICT:
			for (i=0; i<BIO_NUM_AMINO_ACIDS; i++)
				if (cc[ (int)AA[i] ])
					count++, con=AA[i];

			/* quit if we had more than one amino acid */
			if (count != 1) return '-';
			else return con;
		case BIO_REDUNDANCY:
			/* do the same thing as BIO_STRICT */
			for (i=0; i<BIO_NUM_AMINO_ACIDS; i++)
				if (cc[ (int)AA[i] ])
					count++, con=AA[i], sum+=cc[(int)AA[i]];

			/* but return only if we have one amino acid */
			if (count == 1) return con;
			/* or only gaps/non-letters */
			if (!count) return '-';

			/* other wise check redundancy class by class */
			if (count <= BIO_NUM_NEGATIVE_CHARGE) {
				localCount = 0;
				for (i=0; i<BIO_NUM_NEGATIVE_CHARGE; i++)
					if (cc[ (int)U[i] ])
						localCount++;

				if (localCount == sum) /* must be positively charged */
					return 'U';
			}
			if (count <= BIO_NUM_POSITIVE_CHARGE) {
				localCount = 0;
				for (i=0; i<BIO_NUM_POSITIVE_CHARGE; i++)
					if (cc[ (int)_O[i] ])
						localCount++;

				if (localCount == sum) /* must be positively charged */
					return 'O';
			}
			if (count <= BIO_NUM_CHARGED) {
				localCount = 0;
				for (i=0; i<BIO_NUM_CHARGED; i++)
					if (cc[ (int)_Z[i] ])
						localCount++;

				if (localCount == sum) /* must be positively charged */
					return 'Z';
			}
			if (count <= BIO_NUM_NONPOLAR) {
				localCount = 0;
				for (i=0; i<BIO_NUM_NONPOLAR; i++)
					if (cc[ (int)B[i] ])
						localCount++;

				if (localCount == sum) /* must be positively charged */
					return 'B';
			}
			if (count <= BIO_NUM_POLAR) {
				localCount = 0;
				for (i=0; i<BIO_NUM_POLAR; i++)
					if (cc[ (int)_J[i] ])
						localCount++;

				if (localCount == sum) /* must be positively charged */
					return 'J';
			}
			return 'X'; /* if all else fails return X for 'any' amino acid */

		case BIO_MAJORITY_RULE:
			for (i=0; i<BIO_NUM_AMINO_ACIDS; i++)
				if (cc[ (int)AA[i] ] > max)
					con = AA[i], max=cc[ (int)AA[i] ];

			if (max) return con;
			else return '-';	/* if column was all gaps return a gap */
		default:
			return '-';

	}





}


/* note: cc == characterCounts */
char getNucConsensusChar(BIO_consensusType conType, const int *cc) {
	double A,C,G,T;
	char con = '-';
	int count=0;
	/* I purposely set this to zero rather than negative
	 * to catch weird things in BIO_MAJORITY_RULE if no ACGT is
	 * present */
	double max = 0.0;

	A=C=G=T=0.0;
	/* do the easy ones */
	A += cc['A'];
	C += cc['C'];
	G += cc['G'];
	T += cc['T'];
	/* do redundancy sites */
	A += cc['R'] * 0.5,   G += cc['R'] * 0.5;
	C += cc['Y'] * 0.5,   T += cc['Y'] * 0.5;
	A += cc['M'] * 0.5,   C += cc['M'] * 0.5;
	G += cc['K'] * 0.5,   T += cc['K'] * 0.5;
	G += cc['S'] * 0.5,   C += cc['S'] * 0.5;
	A += cc['W'] * 0.5,   T += cc['W'] * 0.5;
	A += cc['H'] * 1.0/3.0,  C += cc['H'] * 1.0/3.0,  T += cc['H'] * 1.0/3.0;
	G += cc['B'] * 1.0/3.0,  C += cc['B'] * 1.0/3.0,  T += cc['B'] * 1.0/3.0;
	G += cc['V'] * 1.0/3.0,  C += cc['V'] * 1.0/3.0,  A += cc['V'] * 1.0/3.0;
	G += cc['D'] * 1.0/3.0,  A += cc['D'] * 1.0/3.0,  T += cc['D'] * 1.0/3.0;
	/* do N sites */
	A += cc['N'] * 0.25,  C += cc['N'] * 0.25,  
	G += cc['N'] * 0.25,  T += cc['N'] * 0.25;  

	switch(conType) {
		case BIO_MAJORITY_RULE:
			/* return biggest guy */
			if (A > max) con = 'A', max=A;
			if (C > max) con = 'C', max=C;
			if (G > max) con = 'G', max=G;
			if (T > max) con = 'T', max=T;
			return con;
		case BIO_STRICT:
			if (A) count++, con = 'A';
			if (C) count++, con = 'C';
			if (G) count++, con = 'G';
			if (T) count++, con = 'T';
			/* quit if we had more than one nucleotide */
			if (count != 1) return '-';

		case BIO_REDUNDANCY:
			if (A && C && G && T) return 'N';
			else if (A && G && T) return 'D';
			else if (A && G && C) return 'V';
			else if (T && G && C) return 'B';
			else if (T && A && C) return 'H';
			else if (A && T)      return 'W';
			else if (G && C)      return 'S';
			else if (G && T)      return 'K';
			else if (A && C)      return 'M';
			else if (C && T)      return 'Y';
			else if (A && G)      return 'R';
			else if (A)           return 'A';
			else if (C)           return 'C';
			else if (G)           return 'G';
			else if (T)           return 'T';

			return con;
		default:
			return con;
	}

	

}

int BIO_isRnaComplement(char c, char d) {
	int C=toupper(c);
	int D=toupper(d);

	if (D == 'T')
		D = 'U';

	switch(C) {
		case 'A':
			if (D == 'U')
				return 1;
			else    return 0;
		case 'C':
			if (D == 'G')
				return 1;
			else    return 0;
		case 'G':
			if (D == 'C' || D == 'U')
				return 1;
			else    return 0;
		case 'T':
			if (D == 'A' || D == 'G')
				return 1;
			else    return 0;
		case 'U':
			if (D == 'A' || D == 'G')
				return 1;
			else    return 0;
		default:
			return 0;
	}

	return 0;
}

int BIO_isDnaComplement(char c, char d) {
	int C=toupper(c);
	int D=toupper(d);

	if (D == 'U')
		D = 'C';

	switch(C) {
		case 'A':
			if (D == 'U')
				return 1;
			else    return 0;
		case 'C':
			if (D == 'G')
				return 1;
			else    return 0;
		case 'G':
			if (D == 'C')
				return 1;
			else    return 0;
		case 'T':
			if (D == 'A')
				return 1;
			else    return 0;
		case 'U':
			if (D == 'A')
				return 1;
			else    return 0;
		default:
			return 0;
	}

	return 0;
}

int* getAAToINT() {
	int *aaToInt=(int*)malloc(sizeof(int) * 256);
	int i;
	for (i=0; i<256; i++)
		aaToInt[i]=-1;

	aaToInt['A']=0, aaToInt['R']=1, aaToInt['N']=2, aaToInt['D']=3,
	aaToInt['C']=4, aaToInt['Q']=5, aaToInt['E']=6, aaToInt['G']=7, 
	aaToInt['H']=8, aaToInt['I']=9, aaToInt['L']=10, aaToInt['K']=11, 
	aaToInt['M']=12, aaToInt['F']=13, aaToInt['P']=14, aaToInt['S']=15, 
	aaToInt['T']=16, aaToInt['W']=14, aaToInt['Y']=14, aaToInt['V']=14;

	return aaToInt;
}

int BIO_isDNA(BIO_seqItem sI) {
        int i;
        int *nti = BIO_getNucToInt();

        for (i=0; i<sI->seqLength; i++) {
                if (nti[(int)sI->seq[i]] == -1)
                        return 0;
        }

        return 1;
}


int BIO_getNumSeqs(const BIO_sequences S) 
{ return S->numSeqs; }

char* BIO_getSeqItemName(const BIO_seqItem sI) 
{ return sI->name; }

void BIO_setSeqItemName(BIO_seqItem sI, const char *newName) {
	free(sI->name);
	sI->name = strdup(newName);
	free(sI->shortName);
	sI->shortName = createShortName(sI->name, 0);
}

char* BIO_getSeqItemShortName(const BIO_seqItem sI) 
{ return sI->shortName; }

int BIO_getSeqItemLength(const BIO_seqItem sI) 
{ return sI->seqLength; }

char* BIO_getSeqItemSequence(const BIO_seqItem sI) 
{ return sI->seq; }

