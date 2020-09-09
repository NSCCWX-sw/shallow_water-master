#include <gptl.h>
#include <stdio.h>
void cc()
{
    printf("haha\n");
}
int main ()
{
    int i,ret;
    ret = GPTLinitialize ();               // Initialize GPTL
    for(i=0;i<10;i++){
        ret = GPTLstart ("main");              // Start a manual timer
        cc();
        ret = GPTLstop ("main");               // Stop the manual timer
    }
    ret = GPTLpr_file ("outfile");         // Write output to "outfile"
    return 0;
}
