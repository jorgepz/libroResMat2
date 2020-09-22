clear all, close all, format compact g

% section properties: vector with the cross-section areas
As = [ 0.001 0.01 ]' ;

% material properties: vector with young moduli values
Es = [ 210e9 ] ;

% --- defintions ---
NodsCoord = [ ...
-10 0 ;
 0 10 ;
 0 0  ;
10 0  ] ;

% inode jnode material section
ElemConec = [ ...
1 2 1 1 ;
2 3 1 2 ; ...
2 4 1 1  ...
] ;

fixeddofs = [ 1 2    5 6 7 8 ] ;

NodalLoads = [ 2 10000 0 ];

# Deformed scale factor
scalefactor = 1e3;

# row vector with the dofs related to supports which are replaced by virtual forces
virtualforcessupports = [ 7 ] ;

# row vector with the truss elements which are replaced by virtual forces
virtualforceselements = [  ] ; 

# degree of freedom whose displacement must be determined (leave empty if none)
unkndispdof = 3 ;
% -----------------------------
