
% section properties: vector with the cross-section areas
As = pi * ( 0.1^2 - (0.01-0.005)^2 ) / 4.0 ;

% material properties: vector with young moduli values
Es = [ 210e9 ] ;

NodsCoord = [ ...
 0 4 ;
 0 0 ;
 4 0 ] ;

% inode jnode material section
ElemConec = [ ...
 1 2 1 1 ;
 2 3 1 1 ; ...
 1 3 1 1   ...
 ] ;

fixeddofs = [ 3 4 5 6 ] ;

NodalLoads = [ 1 50e3 0 ];

# Deformed scale factor
scalefactor = 1e3;

# row vector with the dofs related to supports which are replaced by virtual forces
virtualforcessupports = [ 4 ] ;

# row vector with the truss elements which are replaced by virtual forces
virtualforceselements = [  ] ;

# degree of freedom whose displacement must be determined (leave empty if none)
unkndispdof = 1 ;
%~ % -----------------------------
