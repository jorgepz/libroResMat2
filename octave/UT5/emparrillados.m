% ----------------------------------------------------------------------
% ejemplo de codigo para analisis de emparrillados usando a metodo matricial
% Ejecucion en GNU-Octave - Octubre 2019 - Jorge Perez Zerpa, Ignacio Suarez, Bruno Bouchard
% ----------------------------------------------------------------------
close all, clear all

% parametros
l  = 2 ;  width = 0.05 ;    E  = 210e9 ; nu = 0.3 ;   P  = 4e3   ;

% modulo de corte
G = E/(2*(1+nu)) ;
% propiedades geometricas de seccion cuadrada
J = 0.141*width^4 ; I  = width^4/12 ; A  = width^2    ;

%         x    z      
Nodes = [ 0   l/2  ; ...
          l   l/2  ; ...
          l     0  ] ;

% conectividad: nodo1 nodo2
Conec = [       1     2 ; ...
                2     3 ] ;

Angles = [ 0 -pi/2 ]' ;
fixeddofs = [ 1 2 3 7 ];

nelems = size( Conec,1);
nnodes = size( Nodes,1);

freedofs  = (1:(3*nnodes));
freedofs(fixeddofs) = [] ;

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

   KL([2 3 5 6], [2 3 5 6]) = E*I * [ 12/(lelem^3)  6/(lelem^2) -12/(lelem^3)  6/(lelem^2) ; ...
                                      6/(lelem^2)  4/(lelem  )  -6/(lelem^2)  2/(lelem  ) ; ...
                                    -12/(lelem^3) -6/(lelem^2)  12/(lelem^3) -6/(lelem^2) ; ...
                                      6/(lelem^2)  2/(lelem  )  -6/(lelem^2)  4/(lelem  ) ] ;

  dofsElem = [ (elemNodes(1)*3-2):(elemNodes(1)*3) (elemNodes(2)*3-2):(elemNodes(2)*3) ] ;

  KG( dofsElem, dofsElem ) += R * KL * R' ;
end

FG      = zeros ( 3*nnodes, 1        ) ;
FG(8) = -P/2 ;

% imposicion de condiciones de contorno
KG(fixeddofs, :) = [] ; KG(:, fixeddofs) = [] ;
FG( fixeddofs  ) = [] ;

% resolucion de sistema
u = KG\FG ;    UG = zeros( 3*nnodes,1);    UG(freedofs ) = u 

%ploteo de parrillado
figure, hold on, grid on
diameterStructure = max(max( Nodes))  - min(min(Nodes)) ;
margen = 0.1 * diameterStructure ;
for i=1:nelems
  nodesElem = Conec(i,:) ;
  p1 = Nodes( nodesElem(1),:) 
  p2 = Nodes( nodesElem(2),:) 
  plot3( [p1(1)  p2(1)], [p1(2)  p2(2)], [0 0] ,'b')
end
view(3)
