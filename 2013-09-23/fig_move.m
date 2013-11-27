function  fig_move(src,event)
%fig_move moves the figure by 50% of the current x-axis view 
%   Detailed explanation goes here

if strcmp(get(src,'SelectionType'),'alt')
    limits = get(gca,'XLim');
    limit_length = limits(2)-limits(1);
    set(gca,'XLim',[(limits(1)+(limits_length/2)) (limits(2)+(limits_length/2))])
elseif str
    limits = get(gca,'XLim');
    set(gca,'XLim',[(limits(1)-(limits_length/2)) (limits(2)-(limits_length/2))])
else 
    disp('you done messed up')

end

