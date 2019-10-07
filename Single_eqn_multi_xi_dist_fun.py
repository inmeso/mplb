"""
Python routine to generate code when same equation is to be used
based on the value of Xi.

Eg. In BGK, velocity ID decides the values of various parameters. Although the 
number of equation is just one and only alpha varies. 
"""

import re
import glob
import User_defined_function as Udf


def GetValVariableUsingRegexPassed(PassedVarName, Text, Regex):

    regexExpr = Regex.replace('VariableName', PassedVarName)
    pattern = re.compile(regexExpr)
    
    NumMatchFound = 0

    for match in pattern.finditer(Text):
        NumMatchFound += 1
        value = match.group(6)

    if NumMatchFound > 1:
        print 'Warning! Multiple value of ',PassedVarName,' found'
    if NumMatchFound == 1:
        return value
    if NumMatchFound <= 0:
        print 'Warning! Cannot find the value of',PassedVarName
#End of function definition GetValueofVariable


def ParseRangeVariable(text):
    
    FirstCurlyBracePos = text.find('{')
    SecondCurlyBracePos = text.find('}', FirstCurlyBracePos)

    RangeVariable = text[FirstCurlyBracePos + 1:SecondCurlyBracePos]
    RangeVariable = RangeVariable.split(':')

    for i in range(0, len(RangeVariable)):
        RangeVariable[i] = RangeVariable[i].strip()
    
    MinRange = RangeVariable[0]
    MaxRange = RangeVariable[1]

    NegationSymbol = text.find('~', SecondCurlyBracePos)
    if (NegationSymbol == -1):
        print 'Please insert the ~ symbol between two pairs of {} and {}'
        print 'Aborting code execution'
        return None 

    ThirdCurlyBracePos = text.find('{', SecondCurlyBracePos)
    FourthCurlyBracePos = text.find('}', ThirdCurlyBracePos)

    ListValNotIncluded = text[ThirdCurlyBracePos + 1:FourthCurlyBracePos]
    ListValNotIncluded = ListValNotIncluded.split('|')
    
    for i in range(0, len(ListValNotIncluded)):
        ListValNotIncluded[i] =ListValNotIncluded[i].strip()
    
    # If there is no element to be excluded, then list is given value as 'None' 
    # and this is used later to generate code. This will also act as Safety check. 
    if (ListValNotIncluded == ['']):
        ListValNotIncluded = 'None'

    return MinRange, MaxRange, ListValNotIncluded





text = 'spacedim = 3;'
userRegex = r'( |\t|\n)*(\bVariableName\b)( |\t|\n)*(=)( |\t|\n)*([\w\.]+)(;)'
value = GetValVariableUsingRegexPassed('spacedim', text, userRegex)

text = r'{1:10}~{2|3|4}'
#text = r'{1:10}~{ }'

# Check to ensure that function didn't return None.
if(ParseRangeVariable(text)):
    MinRange, MaxRange, ListValNotIncluded = ParseRangeVariable(text)
    print MinRange, MaxRange, ListValNotIncluded

# text  = Udf.InsertForLoop('i', '1', '10', text)
