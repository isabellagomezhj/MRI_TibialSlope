function [] = plotcs_colour(T,length,c)

% Plots a coordinate system defined by a 4x4 local to global tranformation
% matrix, T.
% length scales the length of the axes.

% ALC 2016
% MW 2020 - added colour

% If length is not specified, set to 1.
if nargin==1
  length = 1;
end

if size(c,1)==1
    c = repmat(c,1,3);
end

hold on
quiver3(T(1,4),T(2,4),T(3,4),length*T(1,1),length*T(2,1),length*T(3,1),'color',c(1,:),'linewidth',2) 
quiver3(T(1,4),T(2,4),T(3,4),length*T(1,2),length*T(2,2),length*T(3,2),'color',c(2,:),'linewidth',2)
quiver3(T(1,4),T(2,4),T(3,4),length*T(1,3),length*T(2,3),length*T(3,3),'color',c(3,:),'linewidth',2)


end