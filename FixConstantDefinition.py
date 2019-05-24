# Copyright 2017 the MPLB team. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.
"""
    @brief   Fix the constant definition
    @author  Jianping Meng
    @details Help the OPS library to recognize the constant definition
    usage: python FixConstantDefinition.py variableName newValue
    The current OPS python translator assumes some function arguments as literal numbers,i.e., not a variable. To fix this, this script can be used before invoking OPS to change all constant variables defined in the "h" files to to its actual value in "CPP" files.
"""
from __future__ import print_function
import glob
import re
import sys


def ReadFiles(fileNames):
    texts = []
    for name in fileNames:
        file = open(name, 'r')
        texts.append(file.read())
        file.close()
    return texts


def GetValueofConstant(variableName, texts):
    regexExpr = r'.*(const)( |\t|\n)+(\bint\b)( |\t|\n)+(\b'
    regexExpr += variableName
    regexExpr += r'\b)( |\t|\n)*(=)( |\t|\n)*([0-9]+)( |\t|\n)*(;)'
    pattern = re.compile(regexExpr)
    numMatched = 0
    for text in texts:
        for match in pattern.finditer(text):
            numMatched += 1
            value = match.group(9)
    if numMatched > 1:
        print('Warning: find multipile definitions for ',variableName)
    if numMatched == 1:
        return value
    if numMatched == 0:
        sys.exit('Warning: cannot find the definition for ' + variableName)


def WriteToFile(text, fileName):
    file = open(fileName, 'w')
    file.write(text)
    file.close()


def ReplaceVariable(variableName, value, text):
    regexExpr = r'(\bops_par_loop\b)'
    pattern = re.compile(regexExpr)
    startPos = []
    endPos = []
    numMatched = 0
    for match in pattern.finditer(text):
        callStart = match.end(0)
        callEnd = text.find(';', callStart)
        indStart = text.find(variableName, callStart, callEnd)
        if (indStart != -1):
            numMatched += 1
            startPos.append(indStart)
            endPos.append(indStart + len(variableName))
    if (numMatched >= 1):
        print("Found ", numMatched, " ops_par_loop calls")
        res = text[:startPos[0]] + value
        for ind in range(1, len(startPos)):
            res += (text[endPos[ind - 1]:startPos[ind]] + value)
        res += text[endPos[-1]:]
        return res
    if numMatched == 0:
        print("Found zero matched")
        return None


if len(sys.argv) > 2:
    VariableName = sys.argv[1]
    Value = sys.argv[2]
    print("We will change all its occurrence of ",VariableName, " to ", Value)
else:
    sys.exit("Please enter a variable name!")
HeadFiles = glob.glob("*.h")
CppFiles = glob.glob("*.cpp")

HeadTexts = ReadFiles(HeadFiles)
CppTexts = ReadFiles(CppFiles)
#Value = GetValueofConstant(VariableName, HeadTexts)
for ind in range(len(CppFiles)):
    print("Working on ", CppFiles[ind], " ...")
    text = ReplaceVariable(VariableName, Value, CppTexts[ind])
    if (text != None):
        WriteToFile(text, CppFiles[ind])
