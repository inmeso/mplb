# Copyright 2017 the MPLB team. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.
# Author: Jianping Meng
# Usage: Automatic translation process

#!/bin/bash
if test -d opsversion
then
    echo "Delete the old files"
    rm -r -f opsversion
fi
mkdir opsversion
cp *.h *.cpp opsversion
cd opsversion
python ../FixConstantDefinition.py SPACEDIM
ops.py lbm2d.cpp boundary.cpp flowfield.cpp model.cpp type.cpp evolution.cpp scheme.cpp
python ../FixKernelDeclarition.py
rm -r -f lbm2d.cpp boundary.cpp flowfield.cpp model.cpp type.cpp evolution.cpp scheme.cpp
cp ../Makefile .
echo "Finished but please add include file manually at this moment" 
