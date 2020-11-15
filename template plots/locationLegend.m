function [ LegLocation ] = locationLegend( y_limits, x_end )
% Perform the position of the legend in order to not overlap the datas

    if abs(y_limits(1) - x_end) < abs(y_limits(2) - x_end)
        LegLocation = 'NorthEast';
    else
        LegLocation = 'SouthEast';
    end


end

