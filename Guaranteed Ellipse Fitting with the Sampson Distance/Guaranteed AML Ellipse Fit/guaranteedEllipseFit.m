%   Function: guaranteedEllipseFit
%
%   This function implements the ellipse fitting algorithm described in
%   Z.Szpak, W. Chojnacki and A. van den Hengel
%   "Guaranteed Ellipse Fitting with the Sampson Distance"
%   Proc. 12th European Conference on Computer Vision. ECCV 
%   Firenze,Italy, oct, 2012
%
%   Parameters:
%
%      t             - an initial seed for parameters [a b c d e f] associated with
%                    the equation   a x^2 + b x y + c y^2 + d x + e y + f = 0
%      data_points   - a 2xN matrix where N is the number of data points	         	
%
%   Returns: 
%
%     a length-6 vector [a b c d e f] representing the parameters of the equation
%    
%     a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%     with the additional result that b^2 - 4 a c < 0.
%
%   See Also: 
%
%    compute_guaranteedellipse_estimates
%    levenbergMarquardtStep
%    lineSearchStep
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 25/7/2012 
function [theta] = guaranteedEllipseFit( t, data_points )

% various variable initialisations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primary loop variable
keep_going = true;
% determines whether we use LineSearchStep or LevenbergMarquadtStep
struct.use_pseudoinverse = false;
% in some case a LevenbergMarquardtStep does not decrease the cost
% function and so the parameters (theta) are not updated
struct.theta_updated = false;
% damping parameter in LevenbergMarquadtStep
struct.lambda = 0.01;
% loop counter (matlab arrays start at index 1, not index 0)
struct.k = 1;
% used to modify the tradeoff between gradient descent and hessian based
% descent in LevenbergMarquadtStep
struct.damping_multiplier = 1.2;
% used in LineSearchStep
struct.gamma = 0.00005;
% number of data points
struct.numberOfPoints = length(data_points);
% barrier term that forces parameters to stay in elliptic region
Fprim = [0 0 2; 0 -1 0; 2 0 0];
struct.F = [Fprim zeros(3,3); zeros(3,3) zeros(3,3)];
struct.I = eye(6,6);
% homotopy parameters that weighs contribution of barrier term
struct.alpha = 1e-3;
% data points that we are going to fit an ellipse to
struct.data_points = data_points;

% various parameters that determine stopping criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum loop iterations 
maxIter = 200;
% step-size tolerance
struct.tolDelta = 1e-7;
% cost tolerance
struct.tolCost = 1e-7;
% parameter tolerance
struct.tolTheta = 1e-7;

% various initial memory allocations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocate space for cost of each iteration
struct.cost = zeros(1,maxIter);
% allocate space for the parameters of each iteration
struct.t = zeros(6,maxIter);
% allocate space for the parameter direction of each iteration
struct.delta = zeros(6,maxIter);
% make parameter vector a unit norm vector for numerical stability
t = t / norm(t);
% store the parameters associated with the first iteration
struct.t(:,struct.k) = t;
% start with some random search direction (here we choose all 1's)
% we can initialise with anything we want, so long as the norm of the
% vector is not smaller than tolDeta. The initial search direction 
% is not used in any way in the algorithm. 
struct.delta(:,struct.k) = ones(6,1);
% main estimation loop
 while (keep_going && struct.k < maxIter)
    
    % allocate space for residuals (+1 to store barrier term) 
    struct.r = zeros(struct.numberOfPoints+1,1);
    % allocate space for the jacobian matrix based on AML component
    struct.jacobian_matrix = zeros(struct.numberOfPoints,6);
    % allocate space for the jacobian matrix based on Barrier component
    struct.jacobian_matrix_barrier = zeros(1,6);
    % grab the current parameter estimates
    t = struct.t(:,struct.k);    
    % residuals computed on data points
    for i = 1:struct.numberOfPoints
        m = data_points(:,i);
        % transformed data point
        ux_j = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
        % derivative of transformed data point
        dux_j =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
        
        % outer product
        A = ux_j * ux_j';
        
        % use identity covs ...
        B = dux_j * dux_j';
        
        tBt = t' * B * t;
        tAt = t' * A * t;
        
        % AML cost for i'th data point
        struct.r(i) = sqrt(tAt/tBt);
        
        % derivative AML component
        M = (A / tBt);
        Xbits = B * ((tAt) / (tBt^2));
        X = M - Xbits;
                                
        % gradient for AML cost function 
        grad = (X*t) / sqrt((tAt/tBt));
        % build up jacobian matrix
        struct.jacobian_matrix(i,:) = grad;        
    end
    
    % Barrier term
    tIt = t' * struct.I * t;
    tFt = t' * struct.F * t;

    
    % add the penalty term 
    struct.r(end) = struct.alpha*(tIt/tFt);
    
    % Derivative barrier component
    N = (struct.I / tFt);
    Ybits = struct.F * ((tIt) / (tFt)^2);
    Y = N-Ybits;
    grad_penalty = 2*struct.alpha*Y*t; 
    struct.jacobian_matrix_barrier(1,:) = grad_penalty;
    
    % Jacobian matrix after combining AML and barrier terms
    struct.jacobian_matrix_full = [struct.jacobian_matrix ; struct.jacobian_matrix_barrier];
    
    % approximate Hessian matrix
    struct.H =  struct.jacobian_matrix_full'* struct.jacobian_matrix_full;
    
    % sum of squares cost for the current iteration
    %struct.cost(k) = 0.5*(struct.r'*struct.r);
    struct.cost(struct.k) = (struct.r'*struct.r);
    
    % If we haven't overshot the barrier term then we use LevenbergMarquadt
    % step
    if (~struct.use_pseudoinverse)
        struct = levenbergMarquardtStep(struct);
    else
        struct = lineSearchStep(struct);
    end
    
    % Check if the latest update overshot the barrier term 
    if (struct.t(:,struct.k+1)' * struct.F * struct.t(:,struct.k+1) <= 0)
        % from now onwards we will only use lineSearchStep to ensure
        % that we do not overshoot the barrier 
        struct.use_pseudoinverse = true;
        struct.lambda = 0;
        struct.t(:,struct.k+1) = struct.t(:,struct.k);
        if (struct.k > 1)
            struct.t(:,struct.k) = struct.t(:,struct.k-1); 
        end
        1;
    % Check for various stopping criteria to end the main loop
    elseif (min(norm(struct.t(:,struct.k+1)-struct.t(:,struct.k)),norm(struct.t(:,struct.k+1)+struct.t(:,struct.k))) < struct.tolTheta && struct.theta_updated)
        keep_going = false;
    elseif (abs(struct.cost(struct.k) - struct.cost(struct.k+1)) < struct.tolCost && struct.theta_updated)
        keep_going = false;
    elseif (norm(struct.delta(:,struct.k+1)) < struct.tolDelta && struct.theta_updated)
        keep_going = false;
    end
   
    struct.k = struct.k + 1;
    
 end

1;
 theta = struct.t(:,struct.k);
 theta = theta / norm(theta);

end

