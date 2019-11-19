#
 # Copyright 2019 United Kingdom Research and Innovation
 #
 # Authors: See AUTHORS
 #
 # Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the #following conditions are met:
 #
 # 1. Redistributions of source code must retain the above copyright notice,
 #    this list of conditions and the following disclaimer.
 # 2. Redistributions in binary form must reproduce the above copyright notice
 #    this list of conditions and the following disclaimer in the documentation
 #    and or other materials provided with the distribution.
 # 3. Neither the name of the copyright holder nor the names of its contributors
 #    may be used to endorse or promote products derived from this software
 #    without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 # ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 # POSSIBILITY OF SUCH DAMAGE.
#

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
