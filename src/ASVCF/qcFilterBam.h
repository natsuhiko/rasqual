#define MAXROW 100
#define MAXSTR 1000
#define MAXCOUNT 10000000
//#define READCOUNTINCREMENT 1000000
#define MAXFIELD 1000
//#define BAM_CIGAR_SHIFT 4
//#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

void usage();
void printRawRead(char *row[MAXROW], int wc, long int qst, long int mst, char qstr, char mstr);
void printSam(char *row[MAXROW], int wc, long int qst, long int mst, char qstr, char mstr);
void parseOutputFormat();
//int passesFilters(char *row[MAXROW], long unsigned int *qst, long unsigned int *mst, char *qstr, char *mstr, struct hash *flagCounts);
void getFieldTagAndValue(char *f, char *tag, char *val);
int findTags(char *row[MAXROW], int wc);
void qcFilterBam(char *inFile);
