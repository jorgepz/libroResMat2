## Auxiliar function - nodes2dofs
# Auxiliar function to perform conversion from node number to correponding degree of freedom 
##

function v= nodes2dofs (u,degree)

v=zeros(degree*length(u),1);

for i=1:length(u)
    v( [ (degree*(i-1) + 1):(degree*i) ] ) = [ (degree*(u(i)-1) + 1):(degree*u(i)) ] ;        
end

