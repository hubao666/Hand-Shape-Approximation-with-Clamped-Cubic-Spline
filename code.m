function [S, n] = compute_spline(x_v, y_v, dy_v, h_v)
    syms x
    n = length(x_v) - 1;
    A = zeros(n+1);
    B = zeros(1, n+1);
    
    % Construct matrix a
    for i = 1:n+1
        for j = 1:n+1
            if i == 1 && j == 1
                A(i, j) = 2 * h_v(i);
            elseif i == n+1 && j == n+1 
                A(i, j) = 2 * h_v(i-1);
            elseif j == i-1 
                A(i, j) = h_v(i-1);
            elseif j == i+1
                A(i, j) = h_v(i);
            elseif i == j
                A(i, j) = 2 * (h_v(i-1) + h_v(i));
            else
                A(i, j) = 0;
            end
        end
    end
    
    % construct vector b
    for i = 1:n+1
        if i == 1
            B(i) = 3 * ((1 / h_v(i)) * (y_v(i+1) - y_v(i)) - dy_v(i));  
        elseif i == n+1 
            B(i) = 3 * (dy_v(end) - (1 / h_v(i-1)) * (y_v(i) - y_v(i-1)));
        else
            B(i) = 3 * ((1 / h_v(i)) * (y_v(i+1) - y_v(i)) - (1 / h_v(i-1)) * (y_v(i) - y_v(i-1)));
        end
    end
    
    % solve for c
    c = A \ B';

    b = zeros(1, n);
    d = zeros(1, n);

    for i = 1:n
        b(i) = (1 / h_v(i)) * (y_v(i+1) - y_v(i)) - (h_v(i) / 3) * (2 * c(i) + c(i+1));
        d(i) = (1 / (3 * h_v(i))) * (c(i+1) - c(i));
    end
    
    % construct s
    S = cell(1, n);
    for i = 1:n
        S{i} = y_v(i) + b(i) * (x - x_v(i)) + c(i) * (x - x_v(i))^2 + d(i) * (x - x_v(i))^3;
    end
end

function slope = forward_diff(y_v, h) 
    slope = (1 / (2*h)) * (-3*y_v(1) + 4*y_v(2) - y_v(3));
end

function slope = backward_diff(y_v, h) 
    slope = (1 / (2*h)) * (y_v(1) - 4*y_v(2) + 3*y_v(3));
end

%%% LEFT hand fingers from left to right 

%% First finger
% x = 1~3
y_v1_left = [23 20 18];
y_v1_right = [7 6 5];
endpoint_slope1_left = forward_diff(y_v1_left, 0.2);
endpoint_slope1_right = forward_diff(y_v1_right, 0.2);
x_v1 = [1  1.2 1.5 2 2.5 3];  
y_v1 = [23 20  17  9 7   5]; 
dy_v1 = [endpoint_slope1_left 0 0 0 0 endpoint_slope1_right]; 
n1 = length(x_v1) - 1;
h_v1 = diff(x_v1); 

% x = 1~5.8
y_v2_left = [23 24 24.5];
y_v2_right = [16 15.2 14.5];
endpoint_slope2_left = forward_diff(y_v2_left, 0.2);
endpoint_slope2_right = forward_diff(y_v2_right, 0.2);
x_v2 = [1  1.5   2    2.5 3    3.5  4    4.5 5  5.5 5.8]; 
y_v2 = [23 24.5  24.6 25  24.8 24.5 23.8 20  18 16  14.5];
dy_v2 = [endpoint_slope2_left 0 0 0 0 0 0 0 0 0 endpoint_slope2_right];
n2 = length(x_v2) - 1;
h_v2 = diff(x_v2);

%% Second finger
% x = 5.1~5.8
y_v3_left = [29 27 25.6];
y_v3_right = [20 18.5 14.5];
endpoint_slope3_left = forward_diff(y_v3_left, 0.2);
endpoint_slope3_right = forward_diff(y_v3_right, 0.2);
x_v3 = [5.1  5.2  5.3  5.4  5.5  5.6    5.67 5.7  5.8];
y_v3 = [29   26   25.6 24.8 21   20.5   20.25  20 14.5   ];
dy_v3 = [endpoint_slope3_left 0 0 0 0 0 0 0 endpoint_slope3_right];
n3 = length(x_v3) - 1;
h_v3 = diff(x_v3);

% x = 5.1~9.7
y_v4_left = [29 29.8 30];
y_v4_right = [20 18 16.8];
endpoint_slope4_left = forward_diff(y_v4_left, 0.2);
endpoint_slope4_right = forward_diff(y_v4_right, 0.2);
x_v4 = [5.1 5.5 5.8  6    7     7.5  8  8.5  8.8  9  9.2  9.3  9.5  9.7];
y_v4 = [29  30  30.2 30.5 30.75 30.6 30 28   26   23 21   18.5 17.5 16.8 ];
dy_v4 = [endpoint_slope4_left 0 0 0 0 0 0 0 0 0 0 0 0 endpoint_slope4_right];
n4 = length(x_v4) - 1;
h_v4 = diff(x_v4);

%% Third finger (middle)
% x = 9.7~10
y_v5_left = [16.8 19 30];
y_v5_right = [16.8 19 27];
endpoint_slope5_left = forward_diff(y_v5_left, 0.2);
endpoint_slope5_right = forward_diff(y_v5_right, 0.2);
x_v5 = [9.7  9.8  9.9  10];
y_v5 = [16.8 18   24   27];
dy_v5 = [endpoint_slope5_left 0 0 endpoint_slope5_right];
n5 = length(x_v5) - 1;
h_v5 = diff(x_v5);

% x = 10~14
y_v6_left = [27 32 32.1];
y_v6_right = [16.15 16.1 16];
endpoint_slope6_left = forward_diff(y_v6_left, 0.2);
endpoint_slope6_right = forward_diff(y_v6_right, 0.2);
x_v6 = [10 10.5 10.6  11   11.5   12   12.3 12.5 12.7  13   13.5   14];
y_v6 = [27 32.3 32.4  32.5 32.4   32.3 32   31.8 31.2  30   18.5   16];
dy_v6 = [endpoint_slope6_left 0 0 0 0 0 0 0 0 endpoint_slope6_right];
n6 = length(x_v6) - 1;
h_v6 = diff(x_v6);


%% Fourth finger
% x = 14~19.7
y_v7_left = [16 16.1 16.45];
y_v7_right = [28.5 27 26];
endpoint_slope7_left = forward_diff(y_v7_left, 0.2);
endpoint_slope7_right = forward_diff(y_v7_right, 0.2);
x_v7 = [14  15  15.5   16  16.5  17    17.5  18 18.8  19   19.2 19.3  19.5  19.7];
y_v7 = [16  18  21     23  28    29.3  29.6  30 29.8  29.6 29.3 29.1  29    26];
dy_v7 = [endpoint_slope7_left 0 0 0 0 0 0 0 0 0 0 0 endpoint_slope7_right];
n7 = length(x_v7) - 1;
h_v7 = diff(x_v7);

% x = 18.2~19.7
y_v8_left = [16 17.5 18.5];
y_v8_right = [23 24 26];
endpoint_slope8_left = forward_diff(y_v8_left, 0.2);
endpoint_slope8_right = forward_diff(y_v8_right, 0.2);
x_v8 = [18.2   18.3   18.4   18.5   18.75   19  19.2   19.25   19.3  19.4   19.5   19.6   19.7];
y_v8 = [16     17     17.5   18.5   19      20  20.5   21      23    23.2   23.8   24.5   26 ];
dy_v8 = [endpoint_slope8_left 0 0 0 0 0 0 0 0 0 0 0 endpoint_slope8_right];
n8 = length(x_v8) - 1;
h_v8 = diff(x_v8);

%% Fifth finger
% x = 18.2~25.8
y_v9_left = [16 12 9.5];
y_v9_right = [17.2 17 16.5];
endpoint_slope9_left = forward_diff(y_v9_left, 0.2);
endpoint_slope9_right = forward_diff(y_v9_right, 0.2);
x_v9 = [18.2  18.25  18.3   18.35  18.45 18.6   19   19.6   20     20.5   20.75   21   21.45   22   22.6   23     24   24.5   25   25.6   25.8];
y_v9 = [16    14     13     12     11.5  11     10   9.85   10.4   11.5   12      13   14      15   16     16.4   17   17.3   17.3 17     16.5];
dy_v9 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
n9 = length(x_v9) - 1;
h_v9 = diff(x_v9);

% x = 20~25.8
y_v10_left = [3 3.15 3.2];
y_v10_right = [14.8 15 16.5];
endpoint_slope10_left = forward_diff(y_v10_left, 0.2);
endpoint_slope10_right = forward_diff(y_v10_right, 0.2);
x_v10 = [20   21   22   22.28  22.5   22.8   23   23.8   24.2   24.6   24.8   25.15   25.55   25.75   25.8];
y_v10 = [3    4    5    6      7      8      8.5  10     11     12     13     14      15      16      16.5];
dy_v10 = [endpoint_slope10_left 0 0 0 0 0 0 0 0 0 0 0 0 0 endpoint_slope10_right];
n10 = length(x_v10) - 1;
h_v10 = diff(x_v10);



[S1, n1] = compute_spline(x_v1, y_v1, dy_v1, h_v1);
[S2, n2] = compute_spline(x_v2, y_v2, dy_v2, h_v2);
[S3, n3] = compute_spline(x_v3, y_v3, dy_v3, h_v3);
[S4, n4] = compute_spline(x_v4, y_v4, dy_v4, h_v4);
[S5, n5] = compute_spline(x_v5, y_v5, dy_v5, h_v5);
[S6, n6] = compute_spline(x_v6, y_v6, dy_v6, h_v6);
[S7, n7] = compute_spline(x_v7, y_v7, dy_v7, h_v7);
[S8, n8] = compute_spline(x_v8, y_v8, dy_v8, h_v8);
[S9, n9] = compute_spline(x_v9, y_v9, dy_v9, h_v9);
[S10, n10] = compute_spline(x_v10, y_v10, dy_v10, h_v10);

% print out the results
disp('S1:');
for i = 1:n1
    fprintf('S1{%d}(x) = %s\n', i, char(S1{i}));
end

disp('S2:');
for i = 1:n2
    fprintf('S2{%d}(x) = %s\n', i, char(S2{i}));
end

disp('S3:');
for i = 1:n3
    fprintf('S3{%d}(x) = %s\n', i, char(S3{i}));
end

disp('S4:');
for i = 1:n4
    fprintf('S4{%d}(x) = %s\n', i, char(S4{i}));
end

disp('S5:');
for i = 1:n5
    fprintf('S5{%d}(x) = %s\n', i, char(S5{i}));
end

disp('S6:');
for i = 1:n6
    fprintf('S6{%d}(x) = %s\n', i, char(S6{i}));
end

disp('S7:');
for i = 1:n7
    fprintf('S7{%d}(x) = %s\n', i, char(S7{i}));
end

disp('S8:');
for i = 1:n8
    fprintf('S8{%d}(x) = %s\n', i, char(S8{i}));
end

disp('S9:');
for i = 1:n9
    fprintf('S9{%d}(x) = %s\n', i, char(S9{i}));
end

disp('S10:');
for i = 1:n10
    fprintf('S10{%d}(x) = %s\n', i, char(S10{i}));
end




colors = lines(10); 

% plot the graph
figure;
hold on;
line_width = 2; 

for i = 1:n1
    fplot(S1{i}, [x_v1(i) x_v1(i+1)], 'Color', colors(1,:), 'LineWidth', line_width);
end
for i = 1:n2
    fplot(S2{i}, [x_v2(i) x_v2(i+1)], 'Color', colors(2,:), 'LineWidth', line_width);
end
for i = 1:n3
    fplot(S3{i}, [x_v3(i) x_v3(i+1)], 'Color', colors(3,:), 'LineWidth', line_width);
end
for i = 1:n4
    fplot(S4{i}, [x_v4(i) x_v4(i+1)], 'Color', colors(4,:), 'LineWidth', line_width);
end
for i = 1:n5
    fplot(S5{i}, [x_v5(i) x_v5(i+1)], 'Color', colors(5,:), 'LineWidth', line_width);
end
for i = 1:n6
    fplot(S6{i}, [x_v6(i) x_v6(i+1)], 'Color', colors(6,:), 'LineWidth', line_width);
end
for i = 1:n7
    fplot(S7{i}, [x_v7(i) x_v7(i+1)], 'Color', colors(7,:), 'LineWidth', line_width);
end
for i = 1:n8
    fplot(S8{i}, [x_v8(i) x_v8(i+1)], 'Color', colors(8,:), 'LineWidth', line_width);
end
for i = 1:n9
    fplot(S9{i}, [x_v9(i) x_v9(i+1)], 'Color', colors(9,:), 'LineWidth', line_width);
end
for i = 1:n10
    fplot(S10{i}, [x_v10(i) x_v10(i+1)], 'Color', colors(10,:), 'LineWidth', line_width);
end

hold off;
xlabel('x');
ylabel('f(x)');
title('Clamped Cubic Spline Plot of My Left Hand');


axis([0 26 0 36]);
axis equal;

grid on;

