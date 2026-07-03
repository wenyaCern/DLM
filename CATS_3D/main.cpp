#include "TApplication.h"
#include "TROOT.h"

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    gROOT->SetBatch(true);
    gROOT->ProcessLine(".x /Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/run_CATS3D.C");
    return 0;
}