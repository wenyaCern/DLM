#include <iostream>

#include "Basics.h"
#include "ExtendedCk.h"

#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"

using namespace std;


int main(int argc, char *argv[]){

    printf("\nHello to this lovely CATS tutorial!\n");
    printf(" To find a bug: continue with this tutorial\n");
    printf("-------------------------------------------\n");

    // RunCorrectLambdaCent30To50Mt2To7(true, true);
    fitDataCoreReso(true, true, 2.0, true);
    // fitData(true, true);
    return 0;
}
