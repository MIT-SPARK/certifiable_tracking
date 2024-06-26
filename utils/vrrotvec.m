function w = vrrotvec(a, b)
% Extract axis, angle rotation difference between two unit vectors.

% make sure unit
a = a/norm(a);
b = b/norm(b);

% axis
axis = cross(a,b)/norm(cross(a,b));

% angle
angle = acos(dot(a,b));


if (sum(isnan(axis)))
    axis = [0,0,1];
end

w = [axis, angle];

end