% =====================================================================
% ---   Example of Finite Element Method for a thin plate ---
% Solves a rectangular plate with vertical normal load
% using cuadrilateral plate elements.
% working on GNU Octave v. 3.6.4
% @ Copyleft, license - Feb-2014 - J.M. Perez Zerpa, P. Castrillo
%
% Instituto de Estructuras y Transporte
% Universidad de la Republica
% Montevideo, Uruguay
% ----------------------------------------------------------------------

clc, clear all, close all,

% ====================   Parameters   =====================


% --- Material parameters ------
E   = 3004160  ; % Young modulus
nu  = 0.0   ;    % Poisson coeficient


% --- Geometrical parameters ---
L   = [ 4 4 ]' ;  % L = [ length in dimension "x", length dimension "y"]
t   = 0.1   ;     % thickness


% --- Mesh parameteres ---
nx  = 20  ;                    % # divisions in direction "x"

ny      = ceil(L(2)/L(1)*nx)   ;  % # divisions in direction "y"

nel     = [ nx ny ]' ;
nnos    = nel + 1          ;  % # of nodes in each direction
nnostot = nnos(1)*nnos(2)  ;  % # of total nodes

% --- Load and boundary conditions -----
q   =  -1 ;  % Distributed load on surface
% Boundary conditions:  0 free 1 displacement fixed 2 displacement and rotation fixed
CondSouth = 1 ;
CondNorth = 1 ;
CondWest  = 0 ;
CondEast  = 0 ;

% --- Nodes coordinates ---
lins1   = linspace( 0 , L(1) , nnos(1) )';
lins2   = linspace( 0 , L(2) , nnos(2) )';

% ====================================================


% ================= Process =================

% --- nodes coordinates matrix ---
coordnodes = [] ;  
for i = 1:nnos(2) 
  coordnodes( ( (nnos(1) * (i-1) + 1) : nnos(1)*i ) ,:) = [ lins1 lins2(i)*ones(nnos(1),1) ] ;
end

% ---- Conectivity matrix ---
Mcon = []; %
for j = 1:nel(2)
  for i = 1:nel(1)
    intri = (i-1)+1+(j-1)*nel(1);
    Mcon( intri , : ) = [(j-1)*nnos(1)+i    (j-1)*nnos(1)+i+1   j*nnos(1)+i+1   j*nnos(1)+i] ;
  end
end

% --- Degrees of freedom ---
NodesSouth  = ( 1:nnos(1) )' ;
NodesEast   = ( nnos(1):nnos(1):nnostot )' ;
NodesNorth  = ( nnos(1)*ny+1:nnostot    )' ;
NodesWest   = ( 1:nnos(1):nnos(1)*ny+1  )' ;

flagfuerzapuntual = 0 ;

% fixed degrees of freedom 
Gradlibfijab=[] ;
%
for i=1:length(NodesSouth)
  switch CondSouth
  case 1
    Gradlibfijab = [ Gradlibfijab  3*NodesSouth(i)-2  3*NodesSouth(i)-1                  ] ;
    %~ Gradlibfijab = [ Gradlibfijab  3*NodesSouth(i)-2   ] ;
  case 2
    Gradlibfijab = [ Gradlibfijab  3*NodesSouth(i)-2  3*NodesSouth(i)-1  3*NodesSouth(i) ] ;
  end
end

Gradlibfijder=[] ;
%
for i=1:length(NodesEast)
  switch CondEast
  case 1
    Gradlibfijder = [ Gradlibfijder  3*NodesEast(i)-2                   3*NodesEast(i) ];
  case 2
    Gradlibfijder = [ Gradlibfijder  3*NodesEast(i)-2  3*NodesEast(i)-1 3*NodesEast(i) ] ;
  end
end

Gradlibfijarr=[] ;
%
for i=1:length(NodesNorth)
  switch CondNorth
  case 1
    Gradlibfijarr = [ Gradlibfijarr  3*NodesNorth(i)-2  3*NodesNorth(i)-1  ] ;
    %~ Gradlibfijarr = [ Gradlibfijarr  3*NodesNorth(i)-2    ] ;
  case 2
    Gradlibfijarr = [ Gradlibfijarr  3*NodesNorth(i)-2  3*NodesNorth(i)-1  3*NodesNorth(i) ];
  end
end

Gradlibfijizq=[] ;
%
for i=1:length(NodesWest)
  switch CondWest
  case 1
    Gradlibfijizq=[Gradlibfijizq 3*NodesWest(i)-2                    3*NodesWest(i)  ] ;
  case 2
    Gradlibfijizq=[Gradlibfijizq 3*NodesWest(i)-2  3*NodesWest(i)-1  3*NodesWest(i) ] ;
  end
end

% fixed degrees of freedom
Gradlibfij = [Gradlibfijab   Gradlibfijder   Gradlibfijarr Gradlibfijizq];
Gradlibfij = unique(Gradlibfij) ;

Gradlibfij

%~ stop


% free degrees of freedom
Gradliblib = 1:3*nnostot ;

% center of the plate node
nodomedio = floor(nnostot/2) +1;
% center of the plate dof
gdlwmedio = nodomedio * 3 - 2 ; 

Gradliblib (Gradlibfij) = [];

% --- Build Elements coordinates vectors
neltotal = nel(1)*nel(2); % total # of elements
for i = 1:neltotal
    Xel(:,i) = (coordnodes( Mcon(i,:) , 1 )) ; 
    Yel(:,i) = (coordnodes( Mcon(i,:) , 2 )) ; 
end

% ------------------   Stiffness matrix   ---------------------

D  = E * t^3 / ( 12 * (1-nu^2) )
KG = sparse( 3*nnostot,3*nnostot);

a  = L(1) / ( nel(1) * 2 ) ;
b  = L(2) / ( nel(2) * 2 ) ;

Ke1 = b / ( 6 * a^3 )  * ...
  [ 6         6*a     0    -6      6*a       0    -3       3*a     0     3      3*a       0  ;
    6*a       8*a^2   0   -6*a     4*a^2     0   -3*a     2*a^2    0    3*a     4*a^2     0  ;
    0          0      0     0       0        0     0        0      0     0       0        0  ;
    -6       -6*a     0     6     -6*a       0     3      -3*a     0    -3     -3*a       0  ;
    6*a      4*a^2    0   -6*a     8*a^2     0   -3*a     4*a^2    0    3*a     2*a^2     0  ;
    0         0       0     0       0        0     0        0      0     0       0        0  ;
    -3      -3*a      0     3     -3*a       0     6     -6*a      0    -6     -6*a       0  ;
    3*a      2*a^2    0   -3*a     4*a^2     0   -6*a     8*a^2    0    6*a    4*a^2      0  ;
    0         0       0     0       0        0     0        0      0     0       0        0  ;
    3        3*a      0    -3      3*a       0    -6       6*a     0     6      6*a       0  ;
    3*a     4*a^2     0   -3*a     2*a^2     0   -6*a      4*a^2   0    6*a     8*a^2     0  ;
    0         0       0     0       0        0     0        0      0     0       0        0 ] ;

Ke2 = a / (6*b^3) * ...
  [ 6       0      6*b       3    0    3*b        -3       0      3*b     -6      0     6*b;
    0       0       0        0    0    0           0       0        0      0      0         0;
	  6*b     0      8*b^2     3*b  0    4*b^2      -3*b     0    2*b^2     -6*b    0     4*b^2;
 	  3       0      3*b       6    0    6*b       -6        0       6*b     -3     0      3*b;
	  0       0       0        0    0     0         0        0          0     0     0       0;
	  3*b     0      4*b^2    6*b   0     8*b^2    -6*b      0       4*b^2   -3*b   0     2*b^2;
	  -3      0     -3*b      -6    0    -6*b       6        0      -6*b     3      0     -3*b;
	  0       0       0        0    0    0          0        0       0       0      0       0;
	  3*b     0      2*b^2     6*b  0    4*b^2    -6*b       0       8*b^2   -3*b   0      4*b^2;
	  -6      0     -6*b      -3    0    -3*b       3        0       -3*b    6      0     -6*b;
	  0       0       0        0    0     0         0        0       0       0      0        0;
	  6*b     0      4*b^2    3*b   0     2*b^2    -3*b      0      4*b^2    -6*b   0     8*b^2 ];

Ke3 = nu / (2*a*b) * ...
  [ 1         a              b              -1     0            -b           1      0          0         -1       -a     0;
     a         0            2*a*b      0      0            0             0      0         0        -a        0      0;
     b         2*a*b       0            -b     0             0           0      0         0          0        0      0;
    -1       0              -b           1     -a            b           -1      a         0         1       0      0;
     0         0                0         -a     0          -2*a*b    a       0          0         0       0      0;
    -b        0               0           b    -2*a*b    0            0      0          0          0       0       0;
 	   1          0               0          -1     a             0            1     -a        -b        -1      0      b;
	   0          0               0           a       0             0          -a     0         2*a*b   0       0       0;
	   0           0              0           0      0            0           -b    2*a*b  0           b       0       0;
	  -1        -a              0           1      0            0            -1    0          b          1       a      -b;
	  -a        0               0           0      0             0            0      0          0         a       0        -2*a*b;
	   0         0              0           0      0            0            b      0          0        -b     -2*a*b    0];


Ke4 = (1-nu)/ (30*a*b) * ...
  [ 21       3*a        3*b     -21        3*a          -3*b         21          -3*a         -3*b        -21        -3*a          3*b;
    3*a      8*a^2      0       -3*a      -2*a^2       0            3*a         2*a^2        0           -3*a      -8*a^2       0;
    3*b      0          8*b^2   -3*b       0          -8*b^2     3*b          0	         2*b^2    -3*b       0             -2*b^2;
    -21      -3*a      -3*b	      21      -3*a          3*b         -21          3*a           3*b         21         3*a           -3*b;
    3*a      -2*a^2     0        -3*a      8*a^2       0             3*a       -8*a^2         0          -3*a       2*a^2       0;
    -3*b     0         -8*b^2     3*b      0              8*b^2    -3*b          0            -2*b^2    3*b       0               2*b^2;
    21       3*a       3*b       -21       3*a        -3*b          21           -3*a          -3*b       -21        -3*a            3*b;
    -3*a     2*a^2      0         3*a     -8*a^2        0            -3*a            8*a^2        0           3*a       -2*a^2       0;
    -3*b     0          2*b^2     3*b      0           -2*b^2	     -3*b        0              8*b^2    3*b          0               -8*b^2;
    -21      -3*a      -3*b       21      -3*a         3*b          -21         3*a            3*b        21            3*a           -3*b;
    -3*a    -8*a^2     0          3*a      2*a^2        0            -3*a     -2*a^2         0           3*a          8*a^2          0;
    3*b      0         -2*b^2    -3*b       0          2*b^2       3*b        0            -8*b^2      -3*b             0                 8*b^2 ] ;


Ke = D * (Ke1 + Ke2 + Ke3 + Ke4 ) ;

[v,d] = eig( Ke )    ;
vall = sort(diag(d)) ;
valorespropiosmaschicos = vall(1:3)

for i=1:neltotal
  gdl = [];
  for j=1:4
    gdl = [ gdl  Mcon(i,j)*3-2 Mcon(i,j)*3-1 Mcon(i,j)*3 ] ;
  end
  KG( gdl , gdl ) = KG( gdl , gdl ) + Ke ;
end


% -----------  Nodal forces vector  -------------

qelem = 4*q*a*b* [1/4   a/12   b/12   1/4   -a/12    b/12   1/4   -a/12   -b/12   1/4   a/12   -b/12 ]' ;
qtot  = zeros(3*nnostot,1);
for i=1:neltotal
  gdl = [];
    for j=1:4
      gdl = [ gdl Mcon(i,j)*3-2 Mcon(i,j)*3-1 Mcon(i,j)*3];
    end
  %
  qtot( gdl ) = qtot( gdl ) + qelem;   
end

% --- if the flag is 1 then a load is added at the center point ---
if flagfuerzapuntual == 1 
  gdlwmedio         = nodomedio * 3 - 2 ;
  qtot( gdlwmedio ) = qtot( gdlwmedio )  - 1 ;
end

KGu = KG ( Gradliblib , Gradliblib ) ;

qtotlib = qtot;
qtotlib( Gradlibfij ) = [] ;   

u     = KGu\qtotlib ;

U = zeros( 3 * nnostot , 1 ) ; 

U(Gradliblib) = u ; 

Vz    = U( 1 : 3 : 3*nnostot-2 ) ;
Titax = U( 2 : 3 : 3*nnostot-1 ) ;
Titay = U( 3 : 3 : 3*nnostot   ) ;

flechmaxnum = min(Vz) 
flechmaxteo = 0.322*q/D
errorrelflemax = ( flechmaxnum - flechmaxteo ) / flechmaxteo 

KGr = KG ( Gradlibfij , Gradliblib ) ;

Rfij = KGr * u ;

R = zeros( 3*nnostot , 1 ) + qtot ;
R(Gradlibfij) = Rfij ;    

fuerzaapoyo = R( gdlwmedio ) 

Fz = R( 1 : 3 : 3*nnostot-2 ) ;
Mx = R( 2 : 3 : 3*nnostot-1 ) ;
My = R( 3 : 3 : 3*nnostot   ) ;

% converts the vector displacements to a matrix
MatVz = zeros( nnos(1) , nnos(2) ) ;

for j=1:nnos(2)
  MatVz( 1:nnos(1) , j ) = Vz( ((j-1)*nnos(1)+1):(j*nnos(1)) ) ;
end
% ======================================================================



% ========================== Plotings ==================================

% mesh graph
figure
hold on, grid on
for i=1:neltotal
  plot( coordnodes( Mcon(i,:) , 1) , coordnodes( Mcon(i,:) , 2) )
end
title('mesh graph')

% deformed graph
figure
surf(MatVz)
title('Deformed graph')
xlabel('x'), ylabel('y')

% =========================================================================
