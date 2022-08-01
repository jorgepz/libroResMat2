% Copyright 2022, Jorge M. Perez Zerpa. 
%
% This file is part of FMTS.
%
% FMTS is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 
%
% FMTS is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FMTS.  If not, see <https://www.gnu.org/licenses/>.


% --- compute lengths and inclination of undeformed elements ---
Lengths   =  sqrt ( sum( ( NodsCoord( ElemConec(:,2),:) ...
                         - NodsCoord( ElemConec(:,1),:) ).^2 , 2 ) ) ;

Angles = atan2( ( NodsCoord( ElemConec(:,2),2) - NodsCoord( ElemConec(:,1),2) ) , ...
                ( NodsCoord( ElemConec(:,2),1) - NodsCoord( ElemConec(:,1),1) ) ) ;

nnodes = size( NodsCoord,1);     nelems = size( ElemConec,1);

nfixeddofs          = length(fixeddofs) ;
freedofs            = 1:(2*nnodes);
freedofs(fixeddofs) = [];

KG = sparse( 2*nnodes,2*nnodes ) ;
for i=1:nelems

  l   = Lengths(i) ;    ang = Angles(i);

  nodi = ElemConec(i,1);    nodj = ElemConec(i,2);

  Eele = Es( ElemConec(i,3) ) ;    Aele = As( ElemConec(i,4) ) ;

  elemdofs = nodes2dofs ( [ nodi nodj ]' ,2) ;

  ca = cos(ang);  sa = sin(ang);

  R = [ ca -sa   0   0 ;
        sa  ca   0   0 ;
         0   0  ca -sa ;
         0   0  sa  ca ];

  KL = Eele*Aele/l   * [  1  0 -1  0 ;
                          0  0  0  0 ;
                         -1  0  1  0 ;
                          0  0  0  0 ] ;
  KGelem = R * KL * R'

  KG(elemdofs,elemdofs) =  KG(elemdofs,elemdofs) + KGelem ;
end

FG = zeros(2*nnodes,1);
for i=1:size(NodalLoads,1)
  aux = nodes2dofs ( NodalLoads(i,1), 2 ) ;
  FG( aux ) = FG( aux ) + NodalLoads(i,2:3)' ;
end


%% 3- Process
% Linear system resolution.
%%

K     = KG;
FGred = FG;

K ( fixeddofs , : )    = [];
K ( : , fixeddofs )    = [];
FGred ( fixeddofs    ) = [];

U = K \ FGred ;

UG = zeros(2*nnodes,1);   UG( freedofs) = U ;

UG

%% 4- Posprocess
% Computation of:
% * *ElemStrains*:  vector of strains per element
% * *ElemStresses*: vector of stress per element
% * *NormalForces*: vector of Normal force per element
%%

ElemStrains  = zeros(nelems,1);
ElemStresses = zeros(nelems,1);
NormalForces = zeros(nelems,1);

for i=1:nelems

  l   = Lengths(i) ;
  ang = Angles(i);

  nodi = ElemConec(i,1);
  nodj = ElemConec(i,2);

  Eele = Es( ElemConec(i,3) ) ;
  Aele = As( ElemConec(i,4) ) ;

  elemdofs = nodes2dofs ( [ nodi nodj ]' ,2) ;

  ca = cos(ang);
  sa = sin(ang);

  R = [ ca -sa   0   0 ;
        sa  ca   0   0 ;
         0   0  ca -sa ;
         0   0  sa  ca ];

  KL = ...
  Eele*Aele/l   * [  1  0 -1  0 ;
                     0  0  0  0 ;
                    -1  0  1  0 ;
                     0  0  0  0 ] ;

  Belemloc         = 1/ Lengths(i) * [ -1 0 1 0] ;

  LocalDispl       =  ( R' * UG( elemdofs ) ) ;
  ElemStrains  (i) =  Belemloc * LocalDispl ;
  ElemStresses (i) =  Eele * ElemStrains(i)  ;
  NormalForces (i) =  Aele * ElemStresses(i) ;

end

NormalForces

%% 5- Plots
% Plots normal forces diagram and deformed structure.
%%

% line width (LW) and markersize (MS)
LW = 3; MS = 4;
loadscfactor = max(  max(NodsCoord(:,1)) - min(NodsCoord(:,1)) , ...
                     max(NodsCoord(:,2)) - min(NodsCoord(:,2)) ) ...
               / max( abs( FG ) ) * 0.1 ;

UGmat = reshape(UG,[2,nnodes])' ;
NodsCoordDef = NodsCoord + scalefactor*UGmat(:,1:2);

figure, grid on, hold on, axis equal
quiver( NodsCoord(:,1), NodsCoord(:,2), FG(1:2:end)*loadscfactor ...
                                      , FG(2:2:end)*loadscfactor ,0,'c',"filled")
for i=1:nelems
  xselem = NodsCoord( ElemConec(i,1:2) , 1 );
  yselem = NodsCoord( ElemConec(i,1:2) , 2 ) ;

  plot( xselem, yselem, 'k--', 'linewidth', LW*0.75, 'markersize', MS );

  elemdofs = nodes2dofs( ElemConec(i,1:2),3) ;

  xselem = NodsCoordDef ( ElemConec(i,1:2) , 1 ) ;
  yselem = NodsCoordDef ( ElemConec(i,1:2) , 2 ) ;

  plot( xselem, yselem, 'b-o', 'linewidth', LW, 'markersize', MS );
end
title('Deformed: black: reference configuration, blue: deformed configuration, cyan: loads.')
xlabel('x'); ylabel('x'); print('deformed.png','-dpng');

figure, hold on
quiver( NodsCoord(:,1), NodsCoord(:,2) , FG(1:2:end)*loadscfactor ...
  , FG(2:2:end)*loadscfactor ,0,'c',"filled")

for i=1:nelems

  xselem = NodsCoord( ElemConec(i,1:2) , 1 ) ;
  yselem = NodsCoord( ElemConec(i,1:2) , 2 ) ;

  elemdofs = nodes2dofs( ElemConec(i,1:2), 2 ) ;

  if     NormalForces(i) >0, colornormalforce='b';
  elseif NormalForces(i) <0, colornormalforce='r';
  else                             colornormalforce='k'; end

  plot( xselem, yselem, [colornormalforce '-o'], 'linewidth', LW, 'markersize', MS );

  text( sum(xselem)*0.5*(1.05), sum(yselem)*0.5, sprintf( '%8.2e', NormalForces(i) ) ...
    ,'color',colornormalforce, 'fontsize', 14);
end
axis equal, xlabel('x'), ylabel('y')
title('Normal Forces: cyan: external loads, red: compression forces, blue: tension forces.')
print('normalforces.png','-dpng');
