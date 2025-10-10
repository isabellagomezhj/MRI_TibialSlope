%Find Transition Step (Kneedle)
%Bobby Kent
%7/8/2015
%
%Function that finds the transition step of a load-displacement curve based
%on the 'Kneedle' algortithm outlined in the following article:
%
%Satopaa V, Albrecht J, Irwin D, Raghavan B. Finding a "Kneedle" in a
%haystack: detecting knee points in system behavior. ICDCSW. 2011: 166-71.

function [t_step] = find_transition_step_kneedle(x,y)

%Normalize the points of the curve to a unit square so that the algorithm
%functions the same way regardless of the magnitude of the values in the
%underlying data
x_sn = (x - min(x))./(max(x) - min(x));
y_sn = (y - min(y))./(max(y) - min(y));

%Find a set of differences between the x- and y-values, i.e., the set of
%differences between the set of points
y_d = y_sn - x_sn;

%Find the local minimum of the resulting set
[~, t_step] = min(y_d);

%% Plot Feature
%Comment this in if you want it to generate a plot of x and y with the
%transition step labled

%Initialize figure
% figure();
% hold on
% 
% %Plot
% plot(x,y,'LineWidth',3);
% plot(x(t_step),y(t_step),'ro','MarkerSize',8,'LineWidth',3);
% title('Kneedle');

end