"""
    This Pyhton code will read a User written file for Initilaistion and will generate a code 
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
    regexExpr = r''+String
    #print regexExpr
    pattern = re.compile(regexExpr)
    numMatched = 0

    for match in pattern.finditer(Text):
        pos = match.start()
        position.append(pos)
        numMatched = numMatched + 1
  
    if(numMatched >= 1):
        #print 'Multiple ',String,' found in the text at positions',position
        return position

    #elif(numMatched==1):
        #print String,' found in the text at position',position
    #    return position
#End of FindPositionStringText


def Comment(line,Translated_Text):
    if len(line)==0:
        Translated_Text = Translated_Text + '\n'
    else:
        Translated_Text += '//' + line + '\n'
    return Translated_Text
#End of Function Comment.
    

def GetValueofVariable(VariableName, Text):
    regexExpr = r'( |\t|\n)*(\b'
    regexExpr = regexExpr + VariableName
    regexExpr = regexExpr + r'\b)( |\t|\n)*(=)( |\t|\n)*(\w+\.\w+)(;)' 
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


#----------------------------------------------------------
# Routine to check number of arguements suuplied to routines 
# such as CompiVeloIdx, SpaIdx etc. 
#----------------------------------------------------------

def CheckNumArgs(FunctionName, ArgsGiven):

    # To Do:- Modify this for 2D code. 
    # USe spacedim and subtract one number of arguements.
    if FunctionName == 'CompoVeloSpaIdx' and ArgsGiven!=5:
        print 'CompoVeloSpaIdx expects 5 Arguements but supplied =', ArgsGiven 

    if FunctionName == 'CompoVeloIdx' and ArgsGiven!=2:
        print 'CompoVeloIdx expects 2 Arguements but supplied =', ArgsGiven 

    if FunctionName == 'CompoMacroSpaIdx' and ArgsGiven!=5:
        print 'CompoMacroSpaIdx expects 5 Arguements but supplied =', ArgsGiven 

    if FunctionName == 'SpaIdx' and ArgsGiven!=3:
        print 'SpaIdx expects 3 Arguements but supplied =', ArgsGiven 

# End of routine to check number of arguements.


#----------------------------------------------------------
# Fucntion to parse arguements of Functions
# such as CompiVeloIdx, SpaIdx etc.
#----------------------------------------------------------

def ParseArguements(FunName, Args):
    
    ArgDic = {}
    if FunName == 'CompoVeloSpaIdx':
        ArgDic['CompoId'] = Args[0]
        ArgDic['VeloId'] = Args[1]
        ArgDic['RelSpaIdx_X'] = Args[2]
        ArgDic['RelSpaIdx_Y'] = Args[3]
        ArgDic['RelSpaIdx_Z'] = Args[4]

    if FunName == 'CompoMacroSpaIdx':
        ArgDic['CompoId'] = Args[0]
        ArgDic['MacroVarId'] = Args[1]
        ArgDic['RelSpaIdx_X'] = Args[2]
        ArgDic['RelSpaIdx_Y'] = Args[3]
        ArgDic['RelSpaIdx_Z'] = Args[4]
    
    if FunName == 'SpaIndex':
        ArgDic['RelSpaIdx_X'] = Args[0]
        ArgDic['RelSpaIdx_Y'] = Args[1]
        ArgDic['RelSpaIdx_Z'] = Args[2]

    if FunName == 'CompoVeloIdx':
        ArgDic['CompoId'] = Args[0]
        ArgDic['VeloId'] = Args[1]

    return ArgDic
#End of function to parse arguements.


#----------------------------------------------------------
# Routine to get variable name. 
# This information might be needed in the final code. 
#----------------------------------------------------------

def ParseText(Text, Positions, TypeVar): 
    for pos in Positions:
        Item = {}

        #Which type of variable to search for such as 
        #Dist_f, Weights, Micro_Vel_, Macro_Vars.
        VarNameStartPos = pos + len(TypeVar)
        VarNameEndPos = Text.find('[', VarNameStartPos)
        VarName  = Text[VarNameStartPos:VarNameEndPos].strip()
        Item['VarName'] = VarName
        
        #Everything inside a square bracket from which information
        #has to be extracted.
        ArgStartPos = VarNameEndPos+1
        ArgEndPos = Text.find(']', VarNameStartPos)
        Arguement =  Text[ArgStartPos:ArgEndPos].strip()
        Item['Arguement'] = Arguement

        #The following two positions will be used to insert the 
        #genrated code at the write place in the translated text.
        Item['StartPosTextInsert'] = pos
        Item['EndPosTextInsert'] = ArgEndPos+1

        #print Text[pos:ArgEndPos+1]

        #Getting Routine name such as CompoVeloIdx, CompoMacroSpaIdx etc.
        FunStartPos = ArgStartPos
        FunEndPos = Text.find('(', FunStartPos)
        FunName = Text[FunStartPos: FunEndPos].strip()
        Item['Function'] = FunName


        #Getting Function Arguments which will help in generating code in HiLeMMs.
        FunArgsStart = FunEndPos + 1
        FunArgsEnd = Text.find(')',FunArgsStart)
        FunArgs = Text[FunArgsStart:FunArgsEnd].strip()
        Item['FunArgs'] = FunArgs
        #Parsed_Text.append([Item])
        

        FunArgs = FunArgs.split(',')
        NumArgsFun = len(FunArgs)
        CheckNumArgs(FunName, NumArgsFun)
        
        for i in range(NumArgsFun):
            FunArgs[i] = FunArgs[i].strip()

        Item['ParsedArgs'] = ParseArguements(FunName, FunArgs)
        #print Item
        Parsed_Text.append(Item)
        #Parsed_Text.append(ParseArguements(FunName, FunArgs))
        #print Parsed_Text

# End of Routine to parse Text.


#----------------------------------------------
# Function to merge two dictionaries.
#----------------------------------------------
def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

# End of merge fnction.
#-----------------------------------------------


#----------------------------------------------------
# Function to generate code for the coordinates in 
# the user defined function.
#----------------------------------------------------

def GenCodeCoordinates(Parsed_Text):

    for i in range(0,len(Parsed_Text)):
        
        CodeCoord = {}
        FunName = Parsed_Text[i]['Function']
        
        if FunName == 'SpaIndex':
            
            VariableName = Parsed_Text[i]['VarName']
            RelPos_X = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_X']
            RelPos_Y = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Y']
            RelPos_Z = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Z']

            if VariableName == 'X':
                CodeCoord['GenCode'] = 'coordinates[OPS_ACC_MD0(0,' + RelPos_X + ',' + RelPos_Y + ',' + RelPos_Z + ')]'

            elif VariableName == 'Y':
                CodeCoord['GenCode'] = 'coordinates[OPS_ACC_MD0(1,' + RelPos_X + ',' + RelPos_Y + ',' + RelPos_Z + ')]'

            elif VariableName == 'Z':
                CodeCoord['GenCode'] = 'coordinates[OPS_ACC_MD0(2,' + RelPos_X + ',' + RelPos_Y + ',' + RelPos_Z + ')]'

            else:
                print 'Variable Name for Spatial Coordinates could only be X, Y or Z'

        Parsed_Text[i] = merge_two_dicts(Parsed_Text[i], CodeCoord)

# End of function to gen coordinate code.
#----------------------------------------------------


#-------------------------------------------------------------------------------------------------
# Function to get the correct index of macro vars to be 
# used in generating code.
# MacroVarNames :- List of all macro vars parsed from User written Cpp file.
# CompoIdMacroVars :- List of Id's all macro vars spec. which component macro vars belongs to.
# MacroVarSearch :- For Which macro var, we are genrating the index.
# CompoIdSearch :- The component ID of macro var being searched.   
#--------------------------------------------------------------------------------------------------

def GetIndexMacroVarsforCodeGen(MacroVarNames, CompoIdMacroVars, MacroVarSearch, CompIdSearch):
    Index = 0
    Found = 'False'
    for i in range(0, len(MacroVarNames)):
        if MacroVarSearch == MacroVarNames[i] and CompIdSearch == CompoIdMacroVars[i]:
            Found = 'True'
            Index = i
            break
    if Found == 'False':
        print 'Could not find the macroscopic variable ', MacroVarSearch, 'for component ', CompIdSearch
    else:
        return Index
            
#End of function to generate the correct index.
#--------------------------------------------------------



#----------------------------------------------------
# Function to generate code for the Macroscopic 
# Variables in the user defined function.
#----------------------------------------------------

def GenCodeMacroVars(Parsed_Text):
    
    for i in range(0,len(Parsed_Text)):
        
        CodeMacroVars = {}
        FunName = Parsed_Text[i]['Function']
        
        if FunName == 'CompoMacroSpaIdx':
            
            MacroVarName = Parsed_Text[i]['ParsedArgs']['MacroVarId']
            ComponentId = Parsed_Text[i]['ParsedArgs']['CompoId']
            RelPos_X = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_X']
            RelPos_Y = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Y']
            RelPos_Z = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Z']

            Index = GetIndexMacroVarsforCodeGen(MacroVarNames, CompoIdMacroVars, MacroVarName, ComponentId)
            CodeMacroVars['GenCode'] = 'macroVars[OPS_ACC_MD1(' + str(Index) + ',' + RelPos_X + ',' + RelPos_Y + ',' + RelPos_Z + ')]'

        Parsed_Text[i] = merge_two_dicts(Parsed_Text[i], CodeMacroVars)

# End of function to gen code for macro vars.
#----------------------------------------------------



#----------------------------------------------------
# Function to generate code for the Microscopic Dist
# Fun in the user defined function.
#----------------------------------------------------

def GenCodeDistFun(Parsed_Text):

    for i in range(0,len(Parsed_Text)):
        
        CodeDistFun = {}
        FunName = Parsed_Text[i]['Function']
        
        if FunName == 'CompoVeloSpaIdx':
            
            ComponentId = Parsed_Text[i]['ParsedArgs']['CompoId']
            VeloId = Parsed_Text[i]['ParsedArgs']['VeloId']
            RelPos_X = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_X']
            RelPos_Y = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Y']
            RelPos_Z = Parsed_Text[i]['ParsedArgs']['RelSpaIdx_Z']

            ComponentNumber = 'Component' + ComponentId
            
            if int(VeloId) >= 0 and int(VeloId) <= UserVarsCpp[ComponentNumber]['LattSize']:
                Index = UserVarsCpp[ComponentNumber]['XiStart'] + int(VeloId)
                CodeDistFun['GenCode'] = 'f[OPS_ACC_MD2(' + str(Index) + ',' + RelPos_X + ',' + RelPos_Y + ',' + RelPos_Z + ')]'
            else:
                print '*************************************************************************'
                print 'Xi index for component ', ComponentId, ' should be between 0 and ',UserVarsCpp[ComponentNumber]['LattSize']-1
                print 'Cannot generate code for distribution function'
                print '*************************************************************************'

        Parsed_Text[i] = merge_two_dicts(Parsed_Text[i], CodeDistFun)

# End of function to gen code for distribution fun.
#--------------------------------------------------------------


#----------------------------------------------------
# Function to generate code for the Weights
#----------------------------------------------------

def GenCodeWeights(Parsed_Text):
    
    for i in range(0, len(Parsed_Text)):
        
        CodeWeights = {}
        FunName = Parsed_Text[i]['Function']
        VariableName = Parsed_Text[i]['VarName']

        if FunName == 'CompoVeloIdx' and VariableName=='':

            #print 'Running code gen for weights'
            
            ComponentId = Parsed_Text[i]['ParsedArgs']['CompoId']
            VeloId = Parsed_Text[i]['ParsedArgs']['VeloId']

            ComponentNumber = 'Component' + ComponentId
            
            if int(VeloId) >= 0 and int(VeloId) <= UserVarsCpp[ComponentNumber]['LattSize']:
                Index = UserVarsCpp[ComponentNumber]['XiStart'] + int(VeloId)
                CodeWeights['GenCode'] = 'WEIGHTS[' + str(Index) + ']'
            else:
                print '*************************************************************************'
                print 'Xi index for component ', ComponentId, ' should be between 0 and ',UserVarsCpp[ComponentNumber]['LattSize']-1
                print 'Cannot generate code for the Weights'
                print '*************************************************************************'

        #print CodeWeights
        Parsed_Text[i] = merge_two_dicts(Parsed_Text[i], CodeWeights)

# End of function to gen code for Weights.
#--------------------------------------------------------------



#----------------------------------------------------
# Function to generate code for the Microscopic 
# velocity i.e. XI.
#----------------------------------------------------

def GenCodeXi(Parsed_Text):

    for i in range(0, len(Parsed_Text)):
        
        CodeXi = {}
        FunName = Parsed_Text[i]['Function']
        SpaceDim = int(UserVarsCpp['SpaceDim'])
        VariableName = Parsed_Text[i]['VarName']

        #print VariableName==''
        if FunName == 'CompoVeloIdx' and VariableName != '':

            #print FunName, VariableName
            
            ComponentId = Parsed_Text[i]['ParsedArgs']['CompoId']
            VeloId = Parsed_Text[i]['ParsedArgs']['VeloId']
            
            ComponentNumber = 'Component' + ComponentId
            
            if int(VeloId) >= 0 and int(VeloId) <= UserVarsCpp[ComponentNumber]['LattSize']:
                Index = UserVarsCpp[ComponentNumber]['XiStart'] + int(VeloId)

                if VariableName == 'Cx':    
                    CodeXi['GenCode'] = 'XI[' + str(Index * SpaceDim) + ']'
                    
                elif VariableName == 'Cy':    
                    CodeXi['GenCode'] = 'XI[' + str(Index * SpaceDim + 1) + ']'

                elif VariableName == 'Cz':    
                    CodeXi['GenCode'] = 'XI[' + str(Index * SpaceDim + 2) + ']'

                else:
                    print 'Micro Velocity Variable Name should be Cx, Cy or Cz.'

            else:
                print '*************************************************************************'
                print 'Xi index for component ', ComponentId, ' should be between 0 and ',UserVarsCpp[ComponentNumber]['LattSize']-1
                print 'Cannot generate code for the Weights'
                print '*************************************************************************'

        #print 'Hi',CodeXi
        Parsed_Text[i] = merge_two_dicts(Parsed_Text[i], CodeXi)

# End of function to gen code for Xi.
#--------------------------------------------------------------


#-----------------------------------------------------------------------------
# Function to extract values from the user written 
# CPP file (Eg.- lbm3d_hilemms 

# FunName: function name which is constant as defined  by the iinterface.
# ArgNum : Number of arguement whose value is to be found. Note: Argname 
# might change and cannot be used directly. 
#------------------------------------------------------------------------------
def ParseCppFile(FileName, FunName, ArgNum):
    ArgValue = []
    CppText = ReadFile(FileName)[0]

    FunCallStartPos =  FindPositionStringText(FunName, CppText)[0]
   
    FunParanthesisStart = CppText.find('(', FunCallStartPos)

    FunParanthesisEnd = CppText.find(')', FunParanthesisStart)
    FunArgs = CppText[FunParanthesisStart+1:FunParanthesisEnd].strip()

    FunArgs = FunArgs.split(',')
    NumArgsFun = len(FunArgs)
        
    for i in range(NumArgsFun):
        FunArgs[i] = FunArgs[i].strip()

    ArgNameFindValue = FunArgs[ArgNum]
    ArgPosCppFile = FindPositionStringText(ArgNameFindValue, CppText)[0]
    ArgValStart = CppText.find('{', ArgPosCppFile)
    ArgValEnd = CppText.find('}', ArgValStart)
    
    if FunName == 'DefineComponents':
        String = CppText[ArgValStart+1:ArgValEnd].replace('"','')
        String = String.split(',')
        for i in range(0,len(String)):
            String[i] = String[i].strip()
        ArgValue.append(String)

    elif FunName == 'DefineMacroVars':
        String = CppText[ArgValStart+1:ArgValEnd].replace('"','')
        String = String.split(',')
        for i in range(0,len(String)):
            String[i] = String[i].strip()
        ArgValue.append(String)

    else:
        ArgValue.append(CppText[ArgValStart+1:ArgValEnd])
    #ArgValue[ArgNameFindValue] = CppText[ArgValStart+1:ArgValEnd]
    
    return ArgValue

# End of function to parse cpp file.
# ---------------------------------------------------



#------------------------------------------------------------------
# Function to parse information {x}^{m} to pow(x,m).
#------------------------------------------------------------------

def Parse_Base_Exponents(Text):
    regExpr= r'(\{)( )*([a-zA-Z0-9_\-,\[\]\(\) ]+)(\})( )*(\^)( )*(\{)( )*([0-9]+)( )*(\})'
    pattern = re.compile(regExpr)
    StartPosExpr = []
    EndPosExpr = []

    Base = []
    Exponent = []
    
    #Loop to search base and exponent in user defined functions.
    for match in pattern.finditer(Text):

        #Base_Text saves the combination of matched groups in above
        #Re to get the exact base.
        Base_Text = ''
        for i in range(2, 4):
            if match.group(i) == None:
                Base_Text += ''
            else:
                Base_Text += match.group(i)
        
        StartIndex = match.start(0)
        EndIndex = match.end(0)
        #print Text[StartIndex:EndIndex]

        StartPosExpr.append(StartIndex)
        EndPosExpr.append(EndIndex)

        Base.append(Base_Text)
        Exponent.append(match.group(10))

    #Loop to insert pow(base, exponent) at appropriate places.
    Temp = Text[:StartPosExpr[0]]
    Temp += 'pow(' + Base[0] + ',' + Exponent[0] + ')'

    for i in range(1,len(Base)):
        Temp += Text[EndPosExpr[i-1]+1:StartPosExpr[i]]
        Temp += 'pow(' + Base[i] + ',' + Exponent[i] + ')'

    Temp += Text[EndPosExpr[-1]+1:]
    Text = Temp

    return Text
# End of function to parse base and exponents.
#------------------------------------------------------------------



FileName = 'Dist_func_eqn.txt'
Text = ReadFile(FileName)[0]

#Parsed Text is the one where we are collecting information from user written file.
Parsed_Text = []

#This dictionary will hold values of variables like SpaceDim, Number of Components etc.
UserVarsCpp = {}

# Translated text would be the one where we will use
# Parsed text date and genarate the code that has to be inserted.
Translated_Text = ''

# Convert user written code for exponents into C++ format (i.e. insert pow(Variable, m) into code).
Text = Parse_Base_Exponents(Text)

#print 'Starting to Parse information from User written file'
PlaceHolder = ['Dist_', 'Micro_Vel_', 'Weights', 'Macro_Vars', 'Coord_']
#PlaceHolder = ['Coord_']
for String in PlaceHolder:
    Positions = FindPositionStringText(String, Text)
    ParseText(Text, Positions, String)
#print 'File parsing complete'

CppFileName = GetValueofVariable('CppFileName', Text)
UserVarsCpp['SpaceDim'] = ParseCppFile(CppFileName, 'DefineCase', 1)[0]

LatticeNames = ParseCppFile(CppFileName, 'DefineComponents', 2)[0]
#print LatticeNames

#Creating a dictionary of type C1:{LatNam:Val, Lattsize:Val}
#Need this info to generate code.
XiStart = 0
for i in range(0,len(LatticeNames)):
    CompDetails = {}
    CompNum = 'Component' + str(i)
    CompDetails['LattName'] = LatticeNames[i]
    CompDetails['LattSize'] = int(LatticeNames[i][3:])
    CompDetails['XiStart'] = XiStart
    CompDetails['XiEnd'] = XiStart + int(CompDetails['LattSize']) - 1
    XiStart = XiStart + int(CompDetails['LattSize'])
    UserVarsCpp[CompNum] = CompDetails


CompoIdMacroVars = ParseCppFile(CppFileName, 'DefineMacroVars', 3)[0]
NumComponents = len(LatticeNames)

MacroVarNames = ParseCppFile(CppFileName, 'DefineMacroVars', 1)[0]
#print MacroVarNames, CompoIdMacroVars

# Counting the starting and ending position of macroscopic variable of 
# each component. 
MacroVarStartPos = 0
for i in range(0, NumComponents):
    NumMacroVarsComp = 0

    for idx in range(0, len(CompoIdMacroVars)):
        if i == int(CompoIdMacroVars[idx]):
            NumMacroVarsComp += 1

    MacroVarEndPos = MacroVarStartPos + NumMacroVarsComp -1
    CompNum = 'Component' + str(i)
    UserVarsCpp[CompNum]['MacroStartPos'] = MacroVarStartPos
    UserVarsCpp[CompNum]['MacroEndPos'] = MacroVarEndPos
    MacroVarStartPos = MacroVarEndPos + 1


#print UserVarsCpp
GenCodeCoordinates(Parsed_Text)
GenCodeMacroVars(Parsed_Text)
GenCodeDistFun(Parsed_Text)
GenCodeWeights(Parsed_Text)
GenCodeXi(Parsed_Text)

#print len(Parsed_Text)
# for Val in Parsed_Text:
#     print Val

#Sorting the list according to the insertion position.
Parsed_Text_Sorted = sorted(Parsed_Text, key=lambda k: k['StartPosTextInsert'])

# for Val in Parsed_Text_Sorted:
#     print Val


for i in range(0, len(Parsed_Text_Sorted)):
    Translated_Text = Translated_Text + Parsed_Text_Sorted[i]['GenCode']
    Start_Pos_Copy_Text = Parsed_Text_Sorted[i]['EndPosTextInsert']
    
    if i < len(Parsed_Text_Sorted)-1:
        End_Pos_Copy_Text = Parsed_Text_Sorted[i + 1]['StartPosTextInsert']
        Translated_Text += Text[Start_Pos_Copy_Text:End_Pos_Copy_Text]
    else:
        Translated_Text += Text[Start_Pos_Copy_Text :]

#print Translated_Text

FileToWrite = 'UDF_Translated.cpp'
WriteToFile(Translated_Text, FileToWrite)  

"""
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
"""