#include <stdio.h>
#include "melpe.h"
// decode
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

    while ( fread(txbuf, 1, 11, fin) == 11 )
    {
        melpe_s(spbuf, txbuf);
        fwrite(spbuf, 2, 540, fout);
    }

    fclose(fin);
    fclose(fout);

    return 0;
}
