#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int getLenFromCigar(char* c1);

void main(int argc, char** argv){
	printf("%d\n", getLenFromCigar(argv[1]));
	//return getLenFromCigar(argv[1]);
}

int getLenFromCigar(char* c1){
    int itr=0;
    char op[10];
    long len;
    int offset=0;
    char num[100];
    long res=0;
    int nchar;
    while((sscanf(c1+offset, "%ld%[MIDNSH=PX]", &len, op))>0){
        nchar = sprintf(num, "%ld", len);
        offset += nchar+1;
        itr++;
        if(op[0]=='D' || op[0]=='N' || op[0]=='P' || op[0]=='M' || op[0]=='=' || op[0]=='X'){res += len;};
    }
    return (int)res;
}
