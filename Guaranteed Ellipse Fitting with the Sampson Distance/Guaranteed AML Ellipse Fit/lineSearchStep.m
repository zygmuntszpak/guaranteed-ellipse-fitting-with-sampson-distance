%   Function: lineSearchStep
%
%   This function is used in the main loop of guaranteedEllipseFit in the process
%   of minimizing an approximate maximum likelihood cost function of an
%   ellipse fit to data.  It computes an update for the parameters
%   representing the ellipse, using the pseudo-inverse of a Gauss-Newton
%   approximation of the Hessian matrix. It then performs an inexact
%   line search to determine how far to move along the update direction
%   so that the resulting fit is still an ellipse. 
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

function struct = lineSearchStep( struct )

  % extract variables from data structure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  t = struct.t(:,struct.k);
  jacobian_matrix = struct.jacobian_matrix;
  jacobian_matrix_barrier = struct.jacobian_matrix_barrier;
  r = struct.r;
  I = struct.I;
  lambda = struct.lambda;
  delta = struct.delta(struct.k);
  tolDelta = struct.tolDelta;
  damping_multiplier = struct.damping_multiplier;
  F = struct.F;
  I = struct.I;
  current_cost = struct.cost(struct.k);
  data_points = struct.data_points;
  alpha = struct.alpha;
  gamma = struct.gamma;
  numberOfPoints = struct.numberOfPoints;
  
  % jacobian (vector), r(t)'d/dtheta[r(t)]
  jacob = [jacobian_matrix ; jacobian_matrix_barrier]'*r;
  tFt = t' * F * t;
  
  % We solve for the new update direction in a numerically careful manner
  % If we didn't worry about numerical stability then we would compute 
  % the first new search direction like this:
  
  % update = - pinv(H)*jacob;
  
  % But close to the barrier between ellipses and hyperbolas we may
  % experience numerical conditioning problems due to the nature of the
  % barrier term itself.  Hence we perform the calculation in a
  % numerically more stable way with 
  
  Z = [((jacobian_matrix'*jacobian_matrix) + lambda*I) tFt^4*(jacobian_matrix_barrier'*jacobian_matrix_barrier) ; I -(tFt)^4*I];
  zz = -[jacob ; zeros(6,1)];
  update = pinv(Z,1e-20)*zz;
  % drop the nuisance parameter components
  update = update(1:6);

  % there is no repeat...until construct so we use a while-do
  frac = 0.5;
  while true
    % compute potential update    
    t_potential = t + frac*update;
    delta = frac*update;
    % halve the step-size
    frac = frac / 2 ;
    % compute new residuals on data points
    cost = 0;
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
        
        tBt = t_potential' * B * t_potential;
        tAt = t_potential'  * A * t_potential;
               
        % AML cost for i'th data point
        cost = cost +  tAt/tBt ;               
    end
    
    % Barrier term
    tIt = t_potential' * I * t_potential;
    tFt = t_potential' * F * t_potential;
    
    % add the barrier term
    cost = cost + (alpha*(tIt/tFt))^2;
    
    % check to see if cost function was sufficiently decreased, and whether
    % the estimate is still an ellipse. Additonally, if the step size
    % becomes too small we stop. 
    if (t_potential'* F * t_potential > 0 && (cost < (1-frac*gamma)*current_cost)  || norm(delta) < tolDelta)
     break;
    end
  end
  
  struct.theta_update = true;
  struct.t(:,struct.k+1) = t_potential / norm(t_potential);
  struct.delta(:,struct.k+1) = delta;
  struct.cost(struct.k+1) = cost;
 
  % return a data structure with all the updates
  struct;

end

