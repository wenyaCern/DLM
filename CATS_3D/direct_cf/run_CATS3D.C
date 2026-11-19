#include "CATStools.h"
#include "CATS3D.h"
#include "../fit/shared/Basics3D.h"

void run_CATS3D() {
    CATS3D Kitty;
    Basics_PiPiCATS3D_Setup(Kitty, true, false);

    TH3D* hCk3D = Basics_CATS3D_GetCkHisto(Kitty, "hCk3D_plane");
    Save_TH3D_to_ROOT("test_CATS3D_noCoul_highl.root", hCk3D, "hCk3D_plane");
}
