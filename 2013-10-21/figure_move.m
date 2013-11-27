function  figure_move(src,event)
%figure_move moves the figure by 50% of the current x-axis view
%   Detailed explanation goes here
if strcmp(event.Key,'rightarrow')==1
    limits = get(gca,'XLim');
    limit_length = limits(2)-limits(1);
    set(gca,'XLim',[(limits(1)+(limit_length/2)) (limits(2)+(limit_length/2))])
    disp('right')
elseif strcmp(event.Key,'leftarrow')==1
    limits = get(gca,'XLim');
    limit_length = limits(2)-limits(1);
    set(gca,'XLim',[(limits(1)-(limit_length/2)) (limits(2)-(limit_length/2))])
    disp('left')
elseif strcmp(event.Key,'downarrow')==1
    limits = get(gca,'XLim');
    limit_length = limits(2)-limits(1);
    set(gca,'XLim',[(limits(1)+(limit_length/4)) (limits(2)-(limit_length/4))])
    disp('Zoom In')
elseif strcmp(event.Key,'uparrow')==1
    limits = get(gca,'XLim');
    limit_length = limits(2)-limits(1);
    set(gca,'XLim',[(limits(1)-(limit_length/4)) (limits(2)+(limit_length/4))])
    disp('Zoom out')
else
    disp('Try using using the left and right arrow keys to move the figure.  Use the up and down arrows to zoom')
end

