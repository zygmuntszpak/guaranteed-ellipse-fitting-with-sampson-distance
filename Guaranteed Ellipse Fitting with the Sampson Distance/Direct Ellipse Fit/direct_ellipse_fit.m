% This code is an implementation of the following paper
%
% R. Halif and J. Flusser
% Numerically stable direct least squares fitting of ellipses
% Proc. 6th International Conference in Central Europe on Computer Graphics and Visualization. WSCG '98
% Czech Republic,125--132, feb, 1998

function a = direct_ellipse_fit(data)
    x = data(1,:)';
    y = data(2,:)';
    
    D1 = [x .^ 2, x .* y, y .^ 2]; % quadratic part of the design matrix
    D2 = [x, y, ones(size(x))];    % linear part of the design matrix
    S1 = D1' * D1;                 % quadratic part of the scatter matrix
    S2 = D1' * D2;                 % combined part of the scatter matrix
    S3 = D2' * D2;                 % linear part of the scatter matrix
    T = - inv(S3) * S2';           % for getting a2 from a1
    M = S1 + S2 * T;               % reduce scatter matrix
    M = [M(3, :) ./2; - M(2, :); M(1, :) ./2]; % premultiply by inv(C1)
    [evec, evalue] = eig(M);       % solve eigensystem
    cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :) .^ 2; % evaluate a'Ca
    al = evec(:, find(cond > 0));  % eigenvector for min. pos. eigenvalue
    a = [al; T * al];              % ellipse coefficients
    a = a/norm(a);
end
