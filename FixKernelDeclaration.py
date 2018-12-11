# Copyright 2017 the MPLB team. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.
"""
    @brief   Fix the declaration issue
    @author  Jianping Meng
    @details Help the OPS library to find the declaration of kernels
    usage: python FixKernelDeclaration.py
    The framework of this code is slightly more complicated that the expectation of the OPS so that the OPS python translator cannot automatically find the declaration of kernels. This script helps to insert declration into the sources file for kernels.
"""
import glob
import re

def ReadFiles(fileNames):
    texts = []
    for name in fileNames:
        file = open(name, 'r')
        texts.append(file.read())
        file.close()
    return texts


def WriteToFile(text, fileName):
    file = open(fileName, 'w')
    file.write(text)
    file.close()


def FindEnd(text):
    pattern = re.compile(r'(#)( |\t)*endif')
    numMatched = 0
    end = 0
    for match in pattern.finditer(text):
        numMatched += 1
        if (match.start(0) > end):
            end = match.start(0)
    if numMatched > 1:
        print 'Warning: find multipile #endif '
        return end
    if numMatched == 1:
        return end
    if numMatched == 0:
        return len(text) - 1


def FindDeclaration(text):
    callStr = r'ops_par_loop_'
    indStart = text.find(callStr)
    if (indStart != -1):
        indEnd = text.find(r')', indStart)
        declaration = 'void ' + text[indStart:indEnd + 1] + ";"
        name = text[indStart + len(callStr):text.find(r'(', indStart + len(callStr))]
        spaceInd = name.find(r' ')
        if (spaceInd != -1):
            name = name[:spaceInd]
        return [name, declaration]
    else:
        return None

HeadFiles = glob.glob("*.h")
KernelFiles = glob.glob("*kernel*.h")
HeadFiles = list(set(HeadFiles) - set(KernelFiles))
CppFiles = glob.glob("MPI/*.cpp")
print HeadFiles
HeadTexts = ReadFiles(HeadFiles)
CppTexts = ReadFiles(CppFiles)
Names = []
Declarations = []

for ind in range(len(CppFiles)):
    decl = FindDeclaration(CppTexts[ind])
    if (decl is not None):
        Names.append(decl[0])
        Declarations.append(decl[1])
for fileInd in range(len(HeadTexts)):
    tails = ''
    print 'Working on ' + HeadFiles[fileInd] + ' ...'
    text = HeadTexts[fileInd]
    for ind in range(len(Names)):
        called = (HeadTexts[fileInd].find(Names[ind]) != -1)
        if called:
            tails += Declarations[ind] + '\n'
    # print tails
    fileEnd = FindEnd(text)
    text = text[:fileEnd] + tails + text[fileEnd:]
    WriteToFile(text, HeadFiles[fileInd])
