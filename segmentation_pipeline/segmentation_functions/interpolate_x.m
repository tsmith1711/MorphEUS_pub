function x_out = interpolate_x(y_in, x_vals, y_vals)
    if ismember(y_in, y_vals) %if don't need to interpolate
        x_out = x_vals(find(y_vals == y_in, 1));
    else %solve for y = mx + b
        upper_index = find(y_vals == min(y_vals(y_vals >= y_in)), 1);
        lower_index = find(y_vals == max(y_vals(y_vals <= y_in)), 1);
        slope = ( y_vals(upper_index)- y_vals(lower_index) ) / ( x_vals(upper_index)-x_vals(lower_index) );
        b = y_vals(upper_index) - (slope * x_vals(upper_index));
 
        x_out = (y_in - b)/slope;
    end
end