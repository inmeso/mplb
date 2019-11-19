"""
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
"""

"""
    This Python code will read a User written file for Initialisation and will generate a code
    that has to be inserted at appropriate place for the purpose of customisable ingredients.
"""

import glob
import re


def ReadFile(fileName):
    text = []
    file = open(fileName, 'r')
    text.append(file.read())
    file.close()
    return text


def WriteToFile(text, fileName):
    file = open(fileName, 'w')
    file.write(text)
    file.close()


def FindPositionStringText(String, Text):

    position = []
    regexExpr = r''+String+'\\b'
    #print regexExpr
    pattern = re.compile(regexExpr)
    numMatched = 0

    for match in pattern.finditer(Text):
        pos = match.start()
        position.append(pos)
        numMatched = numMatched + 1

    if(numMatched > 1):
        #print 'Multiple ',String,' found in the text at positions',position
        return position

    elif(numMatched==1):
        #print String,' found in the text at position',position
        return position
#End of FindPositionStringText


def Comment(line,Translated_Text):
    if len(line)==0:
        Translated_Text = Translated_Text + '\n'
    else:
        Translated_Text += '//' + line + '\n'
    return Translated_Text
#End of Function Comment.


def GetValueofVariable(VariableName, Text):
    regexExpr = r'( |\t|\n)+(\b'
    regexExpr = regexExpr + VariableName
    regexExpr = regexExpr + r'\b)( |\t|\n)*(=)( |\t|\n)*([0-9]+)(;)'
    pattern = re.compile(regexExpr)
    NumMatchFound = 0

    for match in pattern.finditer(Text):
        NumMatchFound += 1
        value = match.group(6)

    if NumMatchFound > 1:
        print 'Warning! Multiple value of ',VariableName,' found'
    if NumMatchFound == 1:
        return value
    if NumMatchFound <= 0:
        print 'Warning! Cannot find the value of',VariableName
#End of function definition GetValueofVariable


def Code(Line, Translated_Text):
    Translated_Text += Line + '\n'
    return Translated_Text
#End of function code defintion.



def InsertInitCodeHilemms(File, TextToInsert, NumberSpaceDim):

    Text = ReadFile(FileName)

    #StartPosInitFunction: will be used at starting position to search for OPS 3D and insert
    #the parsed text.
    StartPosInitFunction = FindPositionStringText('KerSetInitialMacroVarsHilemms', Text[0])


    #To Do :- We can limit the search to KerSetInitialMacroVarsHilemms as a safeguard
    #against destroying the original code.

    if NumberSpaceDim == '2':
        StartPosition = Text[0].find('#ifdef OPS_2D', StartPosInitFunction[0])

    elif NumberSpaceDim == '3':
        StartPosition = Text[0].find('#ifdef OPS_3D', StartPosInitFunction[0])

    EndPosition = Text[0].find('#endif', StartPosition)
    #print StartPosition, EndPosition

    Temp = Text[0]
    Text[0] = Temp[0:StartPosition+len('#ifdef OPS_3D')]
    Text[0] += TextToInsert
    Text[0] += Temp[EndPosition:]

    WriteToFile(Text[0], File)
#End of routine Insert Initialisation code for Hilemms.


FileName = 'Initialisation.txt'

Text = ReadFile(FileName)
Translated_Text = ''

#Read the Universal constants.
StartPosition = FindPositionStringText('#Universal_constants', Text[0])
EndPosition = FindPositionStringText('#End_Universal_constants', Text[0])
Translated_Text = Comment('', Translated_Text)
Translated_Text = Comment('Universal_constants',Translated_Text)
Translated_Text = Translated_Text + Text[0][StartPosition[0]+len('#Universal_constants'):EndPosition[0]] + '\n'
Translated_Text = Comment('End_Universal_constants',Translated_Text)

#Insert Blank Lines for Readability.
Translated_Text = Comment('',Translated_Text)
Translated_Text = Comment('',Translated_Text)

#Read the User Defined Constants.
StartPosition = FindPositionStringText('#User_defined_constants', Text[0])
EndPosition = FindPositionStringText('#End_User_defined_constants', Text[0])
#print StartPosition, EndPosition
Translated_Text = Comment('User_defined_constants',Translated_Text)
Translated_Text = Translated_Text + Text[0][StartPosition[0]+len('#User_defined_constants'):EndPosition[0]] + '\n'
Translated_Text = Comment('End_User_defined_constants',Translated_Text)

#Read the value of number of spacial dimensions.
NumberSpaceDim = GetValueofVariable('SpaceDim', Text[0])
#print 'The Current case has',NumberSpaceDim,' number of spatial dimensions '


#---------------------------------------------------------------------------------------
#Write the coordinates information that has to be inserted into the HiLeMMS code.
#---------------------------------------------------------------------------------------
Translated_Text = Comment('',Translated_Text)
Translated_Text = Comment('',Translated_Text)

Translated_Text = Code('Real X;',Translated_Text)
Translated_Text = Code('Real Y;',Translated_Text)
if NumberSpaceDim == '3':
    Translated_Text = Code('Real Z;',Translated_Text)

Translated_Text = Code('X = coordinates[OPS_ACC_MD0(0, 0, 0, 0)];',Translated_Text)
Translated_Text = Code('Y = coordinates[OPS_ACC_MD0(1, 0, 0, 0)];',Translated_Text)
if NumberSpaceDim == '3':
    Translated_Text = Code('Z = coordinates[OPS_ACC_MD0(2, 0, 0, 0)];',Translated_Text)
#-----------------------------------------------------------------------------------------


#Replace all Coord_x or Coord_X by X.
Coordinates = ['x', 'X', 'y', 'Y', 'z', 'Z']
for i in range(len(Coordinates)):
    VariabletoReplace = 'Coord_' + Coordinates[i]
    if i<=1:
        Text[0] = Text[0].replace(VariabletoReplace, 'X')
    elif i==2 or i==3:
        Text[0] = Text[0].replace(VariabletoReplace, 'Y')
    else:
        Text[0] = Text[0].replace(VariabletoReplace, 'Z')
#print Text[0]


#Parse Coord_X^{m}---->X^{m} and insert pow(x,m)
regExpr = r'([X-Z])\^{'
pattern = re.compile(regExpr)
StartPosExpr = []
EndPosExpr = []

Base = []
Exponent = []

itr = 0     #Loop iterator.

#Loop to search base and exponent in user defined dunctions.
for match in pattern.finditer(Text[0]):

    StartIndex = match.start(0)
    EndIndex = Text[0].find('}',StartIndex)

    #print Text[0][StartIndex]
    #print Text[0][EndIndex]

    StartPosExpr.append(StartIndex)
    EndPosExpr.append(EndIndex)

    Base.append(Text[0][match.start(0)])
    Exponent.append(Text[0][StartIndex+3:EndIndex])

#print StartPosExpr, EndPosExpr, Base, Exponent

#Loop to insert pow(base, exponent) at appropriate places.
Temp = Text[0][:StartPosExpr[0]]
Temp += 'pow(' + Base[0] + ',' + Exponent[0] + ')'

for i in range(1,len(Base)):
    Temp += Text[0][EndPosExpr[i-1]+1:StartPosExpr[i]]
    Temp += 'pow(' + Base[i] + ',' + Exponent[i] + ')'

Temp += Text[0][EndPosExpr[-1]+1:]
Text[0] = Temp
#print Temp


#-------------------------------------------------------------------------------
#Insert the macroscopic variable code on LHS of user written expression.
#Eg. Macro_rho will be replaced by macroVars[OPS_ACC_MD2(0, 0, 0, 0)].
#-------------------------------------------------------------------------------

regExpr = r'Macro_'
pattern = re.compile(regExpr)
StartPosExpr = []
EndPosExpr = []

for match in pattern.finditer(Text[0]):

    StartIndex = match.start(0)
    EndIndex = Text[0].find('=',StartIndex)

    StartPosExpr.append(StartIndex)
    EndPosExpr.append(EndIndex)

Temp = Text[0][:StartPosExpr[0]] + 'macroVars[OPS_ACC_MD2(0, 0, 0, 0)]'
for i in range(1,len(StartPosExpr)):
    Temp += Text[0][EndPosExpr[i-1]:StartPosExpr[i]]
    Temp += 'macroVars[OPS_ACC_MD2(' + str(i) + ', 0, 0, 0)]'

Temp += Text[0][EndPosExpr[-1]:]
Text[0] = Temp
#print Text[0]
#-------------------------------------------------------------------------------


#Insert Blank Lines for Readability.
Translated_Text = Comment('',Translated_Text)
Translated_Text = Comment('',Translated_Text)


#Read and Insert the modified function lines.
StartPosition = FindPositionStringText('#Insert the Function/Formula for initialisation here', Text[0])
EndPosition = FindPositionStringText('#End Function for initialistion', Text[0])
Translated_Text = Comment('User Defined Function for Initialisation',Translated_Text)
Translated_Text = Comment('',Translated_Text)
Translated_Text = Translated_Text + Text[0][StartPosition[0]+len('#Insert the Function/Formula for initialisation here.'):EndPosition[0]] + '\n'
Translated_Text = Comment('End Function for initialistion',Translated_Text)

FileToWrite = 'Translated_Initialisation.cpp'
WriteToFile(Translated_Text, FileToWrite)  #Text is a list where first value is the entite text to be written.


FileName = 'hilemms_ops_kernel.h'
InsertInitCodeHilemms(FileName, Translated_Text, NumberSpaceDim)