function DetrapHW4(Y,f,x0,y0,x_end,epsilon,h,h_min,h_max)
%%% Homework # 3, 2019310290 Sangman Jung

% 1. Remark: The problem being solved is Y'=f(x,Y), Y(x0)=y0, for x0 =< x
% =<x_end, using the method described earlier in the section. The
% approximate solution values are printed at each node point. The error
% parameter epsilon and the stepsize parameters were discussed earlier in
% the section. The variable ier is an error indicator, output when exiting
% the algorithm: ier = 0 means a normal return; ier = 1 means that the
% integration was terminated due to a necessary h < h_min.

% Table 6.9 Example of algorithm Detrap
fprintf('Table 6.9 Example of algorithm Detrap\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('|   x_{n}  |     h    |  y_{n}     |  Y(x_{n})-y_{n}  |    trunc   |\n');
fprintf('--------------------------------------------------------------------\n');

% 2. Initialize:
loop = 1; ier = 0;

% 3. Remark: Choose an initial value of h
while 1
    % 4. Calculate y_{h}(x0+h), y_{h/2}(x0+(h/2)), y_{h/2}(x0+h) using
    % method (6.5.2). In each case, use the Euler predictor (6.5.12) and
    % follow it by two iterations of (6.5.3).
    y_p = y0 + h*f(x0,y0);
    y_p = y0 + (h/2)*(f(x0,y0)+f(x0+h,y_p));
    y_h = y0 + (h/2)*(f(x0,y0)+f(x0+h,y_p)); % y_{h}(x0+h)
    y_p = y0 + (h/2)*f(x0,y0);
    y_p = y0 + (h/4)*(f(x0,y0)+f(x0+h/2,y_p));
    y_h22 = y0 + (h/4)*(f(x0,y0)+f(x0+h/2,y_p)); % y_{h/2}(x0+(h/2))
    y_p = y_h22 + (h/2)*f(x0+h/2,y_h22);
    y_p = y_h22 + (h/4)*(f(x0+h/2,y_h22)+f(x0+h,y_p));
    y_h2 = y_h22 + (h/4)*(f(x0+h/2,y_h22)+f(x0+h,y_p)); % y_{h/2}(x0+h)
    
    % 5. For the error in y_{h}(x0+h), use
    trunc = (4/3)*(y_h2-y_h);
    
    % 6. If epsilon*h/4 =< abs(trunc) =< epsilon*h, or if loop = 2, then x1
    % = x0 + h, y1 = y_{h}(x0+h), print x1, y1, and go to step 10.
    if (epsilon*h/4 <= abs(trunc) && abs(trunc) <= epsilon*h) || loop == 2
        x1 = x0 + h;
        y1 = y_h;
        fprintf('|  %1.4f  |  %1.4f  |  %1.6f  |     %1.2e     |  %1.2e  |\n',[x1 h y1 Y(x1)-y1 trunc]);
        break
    end
    
    % 7. Calculate D_3y = Y^{3}(x0) from (6.6.4). If D_3y ~= 0, then
    D_3y = (f(x0+h,y_h2)-2*f(x0+(h/2),y_h22)+f(x0,y0))/(h^2/4);
    if D_3y == 0
        h = h_max;
        loop = 2;
    else
        h = sqrt((6*epsilon)/abs(D_3y));
    end
    
    % 8. If h < h_min, then ier = 2 and exit. If h > h_max, then h = h_max,
    % ier = 1, loop = 2.
    if h < h_min
        ier = 2;
        return
    end
    if h > h_max
        h = h_max;
        ier = 1;
        loop = 2;
    end
    % 9. Go to step 4.
end

% 10. Remark : This portion of the algorithm contains the regular
% predictor-corrctor step with error control.
while 1
    while 1
        % 11. Let x2 = x1 + h, and y2^{0} = y0 + 2*h*f(x1,y1). Iterate
        % (6.5.3) once to obtain y2.
        x2 = x1 + h;
        y_p = y0 + 2*h*f(x1,y1);
        y2 = y1 + (h/2)*(f(x1,y1)+f(x2,y_p));
        
        % 12. trunc = -(y2 - y2^{0})/6
        trunc = -(y2 - y_p)/6;
        
        % 13. If abs(trunc) > epsilon*h or abs(trunc) < epsilon*h/4, then
        % go to step 16.
        if abs(trunc) < ((epsilon*h)/4) || epsilon*h < abs(trunc)
            break
        end
        
        % 14. Print x2, y2.
        if round(x2,4) == 0.2725 || round(x2,4) == 1.9595 || round(x2,4) == 7.6959
            fprintf('--------------------------------------------------------------------\n');
        end
        fprintf('|  %1.4f  |  %1.4f  |  %1.6f  |     %1.2e     |  %1.2e  |\n',[x2 h y2 Y(x1)-y1 trunc]);
        
        % 15. x0 = x1, x1 = x2, y0 = y1, y1 = y2. If x1 < x_end, then go to
        % step 11. Otherwise exit.
        x1 = x2; y0 = y1; y1 = y2;
        if ~(x1 < x_end)
            return
        end
    end
    % 16. Remark : Change the stepsize.
    % 17. x0 = x1, y0 = y1, h0 = h, and calculate h using (6.6.8)
    x0 = x1; y0 = y1; h0 = h;
    h = sqrt((epsilon*(h0^3))/(2*abs(trunc)));
    
    % 18. h = min{h, 2*h0}
    h = min(h,2*h0);
    
    % 19. If h < h_min, then ier = 2 and exit. If h > h_max, then ier = 1
    % and h = h_max.
    if h < h_min
        ier = 2;
        return
    end
    if h > h_max
        ier = 1;
        h = h_max;
    end
    
    % 20. y1^{0} = y0 + h*f(x0,y0), and iterate twice in (6.5.3) to
    % calculate y1. Also, x1 = x0 + h.
    y_p = y0 + h*f(x0,y0);
    y_p = y0 + (h/2)*(f(x0,y0)+f(x0+h,y_p));
    y1 = y0 + (h/2)*(f(x0,y0)+f(x0+h,y_p));
    x1 = x0 + h;
    
    % 21. Print x1, y1.
    if round(x1,4) == 0.2725 || round(x1,4) == 1.9595 || round(x1,4) == 7.6959
        fprintf('--------------------------------------------------------------------\n');
    end
    fprintf('|  %1.4f  |  %1.4f  |  %1.6f  |     %1.2e     |  %1.2e  |\n',[x1 h y1 Y(x1)-y1 trunc]);
    
    % 22. If x1 < x_end, then go to step 10. Otherwise, exit.
    if ~(x1 < x_end)
        return
    end
end