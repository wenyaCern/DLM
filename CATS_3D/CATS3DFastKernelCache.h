#ifndef CATS3D_FAST_KERNEL_CACHE_H
#define CATS3D_FAST_KERNEL_CACHE_H

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS3D.h"

#include "TString.h"

bool SaveCATS3DFastKernel(CATS3D& kitty,
                          const TString& fileName,
                          const TString& treeName="CATS3DFastKernel6D",
                          const bool buildIfNeeded=true);

bool LoadCATS3DFastKernel(CATS3D& kitty,
                          const TString& fileName,
                          const TString& treeName="CATS3DFastKernel6D");

bool SaveCATS3DFastKernelByQ(CATS3D& kitty,
                             const TString& fileName,
                             const TString& treeName="CATS3DFastKernel6DByQ",
                             const bool buildIfNeeded=true);

bool LoadCATS3DFastKernelByQ(CATS3D& kitty,
                             const TString& fileName,
                             const TString& treeName="CATS3DFastKernel6DByQ");

#endif
