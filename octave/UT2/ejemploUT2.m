
% add src folder to path
clear all, close all

problem_name = 'Example_1'

% section properties: vector with the cross-section areas
As = [.001 .01 ]' ;

% material properties: vector with young moduli values
Es = [ 210e9 ] ;

% coordinates (with origin in node B)
NodsCoord = [  0  0  ;
               10 10 ;
               10 0  ;
               20 0  ] ;

% connectiviy: i-node j-node material section
ElemConec = [ 1 2 1 1 ;
              3 2 1 2 ; ...
              4 2 1 1 ] ;

fixeddofs = [ 1 2 5 6 7 8 ] ;

% matrix with nodal info per row:
%             node fx fy
NodalLoads = [ 2 10e3 0 ];

# row vector with the dofs related to supports which are replaced by virtual forces
virtualforcessupports = [  ] ;

# row vector with the truss elements which are replaced by virtual forces
virtualforceselements = [ 3 ] ;

# degree of freedom whose displacement must be determined (leave empty if none)
unkndispdof = 3 ;
% -----------------------------

scalefactor = 1e3;

FEMTrusS
