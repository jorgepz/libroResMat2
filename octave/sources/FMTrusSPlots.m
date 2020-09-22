LW = 2;
MS = 5;
loadscfactor = max(  max(NodsCoord(:,1)) - min(NodsCoord(:,1)) , ...
                     max(NodsCoord(:,2)) - min(NodsCoord(:,2)) ) ...
               / max( abs(Fext) ) * 0.2 ;

figure, hold on
quiver( NodsCoord(:,1), NodsCoord(:,2) , Fext(1:2:end)*loadscfactor ...
  , Fext(2:2:end)*loadscfactor ,0,'g',"filled")
quiver( NodsCoord(:,1), NodsCoord(:,2) , Rext(1:2:end)*loadscfactor ...
  , Rext(2:2:end)*loadscfactor ,0,'r',"filled")

for i=1:nelems

  xselem = NodsCoord( ElemConec(i,1:2) , 1 ) ;
  yselem = NodsCoord( ElemConec(i,1:2) , 2 ) ;
  
  plot( xselem, yselem, 'k--o', 'linewidth', LW, 'markersize', MS );
  
  elemdofs = nodes2dofs( ElemConec(i,1:2), 2 ) ;

  if     ResultNormalForces(i) >0, colornormalforce='b';
  elseif ResultNormalForces(i) <0, colornormalforce='r';
  else                             colornormalforce='k'; end
  
  text( sum(xselem)*0.5, sum(yselem)*0.5, sprintf( '%8.2e', ResultNormalForces(i) ) ...
    ,'color',colornormalforce, 'fontsize', 14);
end
axis equal, xlabel('x'), ylabel('y')
title('Forces results: green: external loads, red: support reactions and blue normal forces.')
