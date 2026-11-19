#include "TApplication.h"
#include "TROOT.h"

#include <filesystem>
#include <string>

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    gROOT->SetBatch(true);
    const auto macro = std::filesystem::path(__FILE__).parent_path() / "run_CATS3D.C";
    gROOT->ProcessLine((std::string(".x ") + macro.string()).c_str());
    return 0;
}
