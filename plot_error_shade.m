function[] = plot_error_shade(hax,x_mean,x_error,x_col,y_mean,y_error,y_col)

x_alpha = 0.2;
y_alpha = 0.1;

if size(x_mean,1)<size(x_mean,2)
    x_mean = x_mean';
end
if size(y_mean,1)<size(y_mean,2)
    y_mean = y_mean';
    y_error = y_error';
end

if ~isempty(y_error)
    x_vector = [x_mean; flipud(x_mean)];
    fill(x_vector, [y_mean+y_error; flipud(y_mean-y_error)],y_col,'edgecolor','none','FaceAlpha',y_alpha,'Parent',hax);
end

if ~isempty(x_error)
    y_vector = [y_mean, flipud(y_mean)];
    fill([x_mean+x_error,flipud(x_mean-x_error)],y_vector,x_col,'edgecolor','none','FaceAlpha',x_alpha,'Parent',hax);
end
end