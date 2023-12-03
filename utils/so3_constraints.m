function c = so3_constraints(R)

% pull out columns
col1 = R(:,1);
col2 = R(:,2);
col3 = R(:,3);

% cols unit length
unitlen = [1.0 - col1'*col1; ...
           1.0 - col2'*col2; ...
           1.0 - col3'*col3];
% cols orthogonal
orth = [col1'*col2;...
        col2'*col3;...
        col3'*col1];
% cols righthanded
righthand = [cross(col1,col2) - col3;...
             cross(col2,col3) - col1;...
             cross(col3,col1) - col2];
% add to matrix!
c = [unitlen; orth; righthand];

end