
% section properties: vector with the cross-section areas
As = 2e-4 ;

% material properties: vector with young moduli values
Es = [ 200e9 ] ;

NodsCoord = [ ...
 0 1 ;
 1/sqrt(3) 0 ;
 1/sqrt(3)+1 1 ;
 1 0 ] ;

% inode jnode material section
ElemConec = [ ...
 1 2 1 1 ;
 1 4 1 1 ; ...
 2 3 1 1 ; ...
 2 4 1 1 ; ...
 4 3 1 1   ...
 ] ;

fixeddofs = [ 1 2 5 6 ] ;

NodalLoads = [ 2 0 -10e3 ; ...
               4 0 -10e3 ];

# Deformed scale factor
scalefactor = 1e3;

# row vector with the dofs related to supports which are replaced by virtual forces
virtualforcessupports = [ ] ;

# row vector with the truss elements which are replaced by virtual forces
virtualforceselements = [ 4 ] ;

# degree of freedom whose displacement must be determined (leave empty if none)
unkndispdof = 4 ;
% -----------------------------
