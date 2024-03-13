px = zeros(9,9);
center = {5,5};
px(center{:}) = -1;

number = 0;
for i = -2:1
    for j = -2:1
        number = number + 1;
        px(center{1}+i, center{2} + j) = number;
    end
end

% Same produces 16 (new) numbers.
% Radius 1: 8
% Radius 2: 25
% Let's drop to a radius 1 (cut new numbers in half, remove dupe)
% 11 is a dupe
% 1,2,3,4,5,9,13,11 should be dropped
keep = [6,7,8,10,12,14,15,16];