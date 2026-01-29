#include <stdio.h>
#include <string.h>
#include "melpe.h"
// encode
int main(int argc, char *argv[])
{

    FILE *fin;
    FILE *fout;

    unsigned char txbuf[11]; //buffer for encoded melpe frame or silency descryptor
    short spbuf[540]; //buffer for accumulate resampled voice up to melpe frame

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s input output\n", argv[0]);
        return -1;
    }

    fin = fopen(argv[1], "r");
    fout = fopen(argv[2], "w");

    if (!fin || !fout)
    {
        fprintf(stderr, "File could not be opened.\n");
        return -1;
    }

    melpe_i();

    int sp_read;
    while ( (sp_read = fread(spbuf, 2, 540, fin)) == 540 )
    {
        melpe_a(txbuf, spbuf);
        fwrite(txbuf, 1, 11, fout);
//        printf("txbuf[10] = %hhu\n", txbuf[10]);
        memset(spbuf, 0, 540 * 2);
    }


    if (sp_read != 540 && sp_read != 0)
    {
        melpe_a(txbuf, spbuf);
        fwrite(txbuf, 1, 11, fout);
    }

    fclose(fin);
    fclose(fout);
}
