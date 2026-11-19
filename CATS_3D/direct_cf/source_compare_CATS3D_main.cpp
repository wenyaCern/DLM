#include "Basics3D.h"

#include "TROOT.h"

int main()
{
    gROOT->SetBatch(true);
    Basics_RunSourceComparison3D("../../inputs/OutputFiles/source_compare_R5.root");
    return 0;
}
