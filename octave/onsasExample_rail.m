% ----------------------------------------------------------------------
% Ejemplo: archivo de modelo numerico de reticulado plano en ONSAS.
% ----------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/../../../onsas/codigo_onsas_repo/' ];
problemName = 'rail' ; addpath( dirOnsas)

E  = 210e9 ;  nu = 0.3  ;
l  = 3     ;  b  = 0.1 ;
P  = 40e4  ;

% MELCS parameters
materialsParams = { [ 0 1 E nu ] } ;
elementsParams  = { [1] ; [2]  } ;
crossSecsParams = { [ 2 b b ]} ;
loadsParams     = { [ 1 1    0 0 0 0 -P 0 ] } ;
springsParams   = { [ inf 0 inf 0 inf 0 ] ; ...
                    [ 0   0 inf 0 inf 0 ] ; ...
                    [ 0   0 inf 0 0   0 ] } ;
Nodes = [   0  0   0 ; ...
            l  0   0 ; ... 
          2*l  0   0 ; ... 
          3*l  0   0 ; ... 
          4*l  0   0 ; ... 
          5*l  0   0 ; ... 
          6*l  0   0 ; ... 
          7*l  0   0 ; ... 
          8*l  0   0 ; ... 
          1*l  0   l ; ... 
          2*l  0   l ; ... 
          3*l  0   l ; ... 
          4*l  0   l ; ... 
          5*l  0   l ; ...
          6*l  0   l ; ...
          7*l  0   l ] ;

Conec = { [ 0 1 0 0 1   1     ] ; ...
          [ 0 1 1 0 3   2     ] ; ...
          [ 0 1 1 0 3   3     ] ; ...
          [ 0 1 1 0 3   4     ] ; ...
          [ 0 1 1 0 3   5     ] ; ...
          [ 0 1 1 0 3   6     ] ; ...
          [ 0 1 1 0 3   7     ] ; ...
          [ 0 1 1 0 3   8     ] ; ...
          [ 0 1 0 0 2   9     ] ; ...
          [ 1 2 0 1 3   1  2  ] ; ... 
          [ 1 2 0 1 3   1  10 ] ; ... 
          [ 1 2 0 1 3   2  10 ] ; ...
          [ 1 2 0 1 3   11 10 ] ; ...
          [ 1 2 0 1 3   2  3  ] ; ...
          [ 1 2 0 1 3   3  11 ] ; ...
          [ 1 2 0 1 3   3  10 ] ; ...
          [ 1 2 0 1 3   3  4  ] ; ...
          [ 1 2 0 1 3   11 4  ] ; ...
          [ 1 2 0 1 3   11 12 ] ; ...
          [ 1 2 0 1 3   4  12 ] ; ...
          [ 1 2 0 1 3   4  13 ] ; ...
          [ 1 2 0 1 3   5  12 ] ; ...
          [ 1 2 0 1 3   13 12 ] ; ...
          [ 1 2 0 1 3   5  4  ] ; ...
          [ 1 2 0 1 3   5  13 ] ; ...
          [ 1 2 0 1 3   6  13 ] ; ...
          [ 1 2 0 1 3   5  14 ] ; ...
          [ 1 2 0 1 3   6  14 ] ; ...
          [ 1 2 0 1 3   6  5  ] ; ...
          [ 1 2 0 1 3   13 14 ] ; ...
          [ 1 2 0 1 3   6  15 ] ; ...
          [ 1 2 0 1 3   14 15 ] ; ...
          [ 1 2 0 1 3   6  7  ] ; ...
          [ 1 2 0 1 3   15 7  ] ; ...
          [ 1 2 0 1 3   16 7  ] ; ...
          [ 1 2 0 1 3   8  7  ] ; ...
          [ 1 2 0 1 3   15 16 ] ; ...
          [ 1 2 0 1 3   8  16 ] ; ...
          [ 1 2 0 1 3   8  9  ] ; ...
          [ 1 2 0 1 3   9  16 ] } ;

plotParamsVector = [ 3 ] ;

ONSAS
