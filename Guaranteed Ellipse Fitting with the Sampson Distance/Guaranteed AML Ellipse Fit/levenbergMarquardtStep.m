%   Function: levenbergMarquardtStep
%
%   This function is used in the main loop of guaranteedEllipseFit in the process
%   of minimizing an approximate maximum likelihood cost function of an
%   ellipse fit to data.  It computes an update for the parameters
%   representing the ellipse, using the method of Levenberg-Marquardt for
%   non-linear optimisation. 
%   See: http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
%
%   Parameters:
%
%      struct     - a data structure containing various parameters
%                   needed for the optimisation process. 	
%
%   Returns: 
%
%     the same data structure 'struct', except that relevant fields have been updated
%
%   See Also: 
%
%    guaranteedEllipseFit
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 25/7/2012 
function struct = levenbergMarquardtStep( struct )

    % extract variables from data structure
    jacobian_matrix = struct.jacobian_matrix;
    jacobian_matrix_barrier = struct.jacobian_matrix_barrier;
    r = struct.r;
    I = struct.I;
    lambda = struct.lambda;
    delta = struct.delta(struct.k);
    damping_multiplier = struct.damping_multiplier;
    F = struct.F;
    I = struct.I;
    t = struct.t(:,struct.k);
    current_cost = struct.cost(struct.k);
    alpha = struct.alpha;
    data_points = struct.data_points;
    numberOfPoints = struct.numberOfPoints;
    
    % compute two potential updates for theta based on different weightings of
    % the identity matrix.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % jacobian (vector), r(t)'d/dtheta[r(t)]
    jacob = [jacobian_matrix ; jacobian_matrix_barrier]'*r;
       
    %tIt = t' * I * t;
    tFt = t' * F * t;
      
    % We solve for the new update direction in a numerically careful manner
    % If we didn't worry about numerical stability then we would compute 
    % the first new search direction like this:
    
    % update_a = - (H+lambda*I)\jacob;
    
    % But close to the barrier between ellipses and hyperbolas we may
    % experience numerical conditioning problems due to the nature of the
    % barrier term itself.  Hence we perform the calculation in a
    % numerically more stable way with 
    
    Z_a = [((jacobian_matrix'*jacobian_matrix) + lambda*I) tFt^4*(jacobian_matrix_barrier'*jacobian_matrix_barrier) ; I -(tFt)^4*I];
    zz_a = -[jacob ; zeros(6,1)];
    update_a = Z_a\zz_a;
    % drop the nuisance parameter components
    update_a = update_a(1:6);
    
    % In a similar fashion, the second potential search direction could be
    % computed like this:
    
    % update_b = - (H+(lambda/v)*I)\jacob
    
    % but instead we computed it with
    Z_b = [((jacobian_matrix'*jacobian_matrix) + (lambda/damping_multiplier)*I) tFt^4*(jacobian_matrix_barrier'*jacobian_matrix_barrier) ; I -(tFt)^4*I];
    zz_b = -[jacob ; zeros(6,1)];
    update_b = Z_b\zz_b;
    % drop the nuisance parameter components
    update_b = update_b(1:6);
    
    % the potential new parameters are then 
    t_potential_a = t + update_a;
    t_potential_b = t + update_b;

    % compute new residuals and costs based on these updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % residuals computed on data points
    cost_a = 0;
    cost_b = 0;
    for i = 1:numberOfPoints
        m = data_points(:,i);
        % transformed data point
        ux_j = [m(1)^2 m(1)*m(2) m(2)^2 m(1) m(2) 1]';
        % derivative of transformed data point
        dux_j =[2*m(1) m(2) 0 1 0 0; 0 m(1) 2*m(2) 0 1 0]';
        
        % outer product
        A = ux_j * ux_j';
        
        % use identity covs ...
        B = dux_j * dux_j';
        
        t_aBt_a = t_potential_a' * B * t_potential_a;
        t_aAt_a = t_potential_a'  * A * t_potential_a;
        
        t_bBt_b = t_potential_b' * B * t_potential_b;
        t_bAt_b = t_potential_b'  * A * t_potential_b;
        
        % AML cost for i'th data point
        cost_a = cost_a +  t_aAt_a/t_aBt_a ;
        cost_b = cost_b +  t_bAt_b/t_bBt_b ;
               
    end
    
    % Barrier term
    t_aIt_a = t_potential_a' * I * t_potential_a;
    t_aFt_a = t_potential_a' * F * t_potential_a;
    
    t_bIt_b = t_potential_b' * I * t_potential_b;
    t_bFt_b = t_potential_b' * F * t_potential_b;

    
    % add the barrier term
    cost_a = cost_a + (alpha*(t_aIt_a/t_aFt_a))^2;
    cost_b = cost_b + (alpha*(t_bIt_b/t_bFt_b))^2;
    
    % determine appropriate damping and if possible select an update
    if (cost_a >= current_cost && cost_b >= current_cost)
        % neither update reduced the cost
        struct.theta_updated = false;
        % no change in the cost
        struct.cost(struct.k+1) = current_cost;
        % no change in parameters
        struct.t(:,struct.k+1) = t;
        % no changes in step direction
        struct.delta(:,struct.k+1) = delta;
        % next iteration add more Identity matrix
        struct.lambda = lambda * damping_multiplier;
    elseif (cost_b < current_cost)
        % update 'b' reduced the cost function
        struct.theta_updated = true;
        % store the new cost
        struct.cost(struct.k+1) = cost_b;
        % choose update 'b'
        struct.t(:,struct.k+1) = t_potential_b / norm(t_potential_b);
        % store the step direction
        struct.delta(:,struct.k+1) = update_b';
        % next iteration add less Identity matrix
        struct.lambda = lambda / damping_multiplier;
    else
        % update 'a' reduced the cost function
        struct.theta_updated = true;
        % store the new cost
        struct.cost(struct.k+1) = cost_a;
        % choose update 'a'
        struct.t(:,struct.k+1) = t_potential_a / norm(t_potential_a);
        % store the step direction
        struct.delta(:,struct.k+1) = update_a';
        % keep the same damping for the next iteration
        struct.lambda = lambda;
    end

    % return a data structure containing all the updates
    struct;
end

