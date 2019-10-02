% ----------------------------------------------------------------------
% ejemplo de codigo para analisis de emparrillados usando a metodo matricial
% Ejecucion en GNU-Octave - Octubre 2019 - Jorge Perez Zerpa
% ----------------------------------------------------------------------

l  = 2     ;
E  = 210e9 ;
nu = 0.3   ;
P  = 4e3   ;

G = 80.7e9
J = 0.141*0.05^4

I  = .05^4/12 ;
A  = .05^2    ;

%         x    z      
Nodes = [ 0   l/2  ; ...
          l   l/2  ; ...
          l     0  ] ;

%        n1 n2 mat sec          
Conec = [ 1 2 ; ...
          2 3 ] ;

Angles = [ 0 -pi/2 ]' ;

fixeddofs = [ 1 2 3 7 ];

freedofs  = (1:(3*nnodes));
freedofs(fixeddofs) = [] ;

nelems = size( Conec,1);
nnodes = size( Nodes,1);

KG      = sparse( 3*nnodes, 3*nnodes ) ;
 
for i = 1:nelems

  alphay = Angles(i);  ca = cos(alphay); sa = sin(alphay);
  
  R = [ ca  0 -sa  0  0   0 ; ...
        0   1   0  0  0   0 ; ...
        sa  0  ca  0  0   0 ; ...
        0   0   0  ca 0 -sa ; ...
        0   0   0  0  1   0 ; ...
        0   0   0  sa 0  ca ] ;
  
  elemNodes = Conec( i,:);
  lelem     = norm( Nodes( elemNodes(2),:) - Nodes( elemNodes(1),:) ) ;
  
  KL = zeros(6,6);
  
  KL([1 4], [1 4])         = G*J/lelem * [ 1 -1 ; ...
                                          -1  1 ] ; 

  KL([2 3 5 6], [2 3 5 6]) = E*I * [ 12/(l^3)  6/(l^2) -12/(l^3)  6/(l^2) ; ...
                                      6/(l^2)  4/(l  )  -6/(l^2)  2/(l  ) ; ...
                                    -12/(l^3) -6/(l^2)  12/(l^3) -6/(l^2) ; ...
                                      6/(l^2)  2/(l  )  -6/(l^2)  4/(l  ) ] ;

  dofsElem = [ (elemNodes(1)*3-2):(elemNodes(1)*3) (elemNodes(2)*3-2):(elemNodes(2)*3) ] ;

  KG( dofsElem, dofsElem ) += R * KL * R' ;
end

FG      = zeros ( 3*nnodes, 1        ) ;
FG(8) = -P/2 ;

% imposicion de condiciones de contorno
KG(fixeddofs, :) = [] ; KG(:, fixeddofs) = [] ;
FG( fixeddofs  ) = [] ;

% resolucion de sistema
u = KG\FG ;
UG = zeros( 3*nnodes,1);
UG(freedofs ) = u 
