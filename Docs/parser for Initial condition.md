This is a sample of the first version of the parser that can be used for Initial condition code.



1. In this parser, I have used some assumptions. I am working with an idea  that we will give a sample format of the file to the user where he can  insert hi code. The sample format has some lines/comments which help to  parse the code.

Consider the following example.

```c++
#Universal_constants
    double Pi = atan(1.0);
#End_Universal_constants
```

These Universal constants are for our use and we can insert as many  constants as we want. Then the user can use them while writing his own  formula.

The lines 

```c++
#Universal_constants
#End_Universal_constants
```

help to search the relevant text in between.

I have also defined an area where a user can insert his own constants and use them in the formula. 



2. I have put a parameter called SpaceDim in the file (which will be  written by the user). This spacedim is used to generate and insert the  relevant coordinates.

For example, when spacedim=3, then only we need to insert code for z coordinate.



3. Initial conditions will rely heavily on the coordinates. I have allowed  the user to use a coordinate in the form Coord_x or Coord_X in the text  and the underlying Python code will then convert it to appropriate  variables used in the HiLemms code.



4. For a coordinate which is raised to some power m, I have used the following rule.

   ​		Coord_X^{m} which gets parsed to pow(X,m) in the generated code.

   

5. For macroscopic variables, user  can write in the form Macro_rho,  Macro_u etc. (We will give a list to choose from macroscopic variables  as is done in the MPLB code). This is then parsed to appropriate  variables in the HiLemms code.



6. I have some ideas on  how to insert the derivatives and multiple components but this will be  different from what you have suggested to me originally. I will explain  this in brief in my next post.

  

I am giving a sample input file and a sample output file for you to check. 



Sample Input

```c++
#Universal_constants
    double Pi = atan(1.0);
#End_Universal_constants


#User_defined_constants
    double U0 = 1.0;
    double U1 = 2.0;
#End_User_defined_constants

#Enter the Number of spacial dimensional here.
    SpaceDim = 3;


#Insert the Function/Formula for initialisation here
    Macro_rho = Coord_x^{2} + 1;
    Macro_U = U0 * sin(2*Pi*Coord_x) * cos(2*Pi*Coord_y) * Coord_z;
    Macro_V = U0 * cos(2*Pi*Coord_x) * Coord_y^{3.5};
    Macro_W = U1 * sin(Coord_x+Coord_y);
#End Function for initialistion
```



Parsed Output

```c++
//Universal_constants

    double Pi = atan(1.0);

//End_Universal_constants


//User_defined_constants

    double U0 = 1.0;
    double U1 = 2.0;

//End_User_defined_constants


Real X;
Real Y;
Real Z;
X = coordinates[OPS_ACC_MD0(0, 0, 0, 0)];
Y = coordinates[OPS_ACC_MD0(1, 0, 0, 0)];
Z = coordinates[OPS_ACC_MD0(2, 0, 0, 0)];


//User Defined Function for Initialisation

    macroVars[OPS_ACC_MD2(0, 0, 0, 0)]= pow(X,2) + 1;
    macroVars[OPS_ACC_MD2(1, 0, 0, 0)]= U0 * sin(2*Pi*X) * cos(2*Pi*Y) * Z;
    macroVars[OPS_ACC_MD2(2, 0, 0, 0)]= U0 * cos(2*Pi*X) * pow(Y,3.5);
    macroVars[OPS_ACC_MD2(3, 0, 0, 0)]= U1 * sin(X+Y);

//End Function for initialistion
```

