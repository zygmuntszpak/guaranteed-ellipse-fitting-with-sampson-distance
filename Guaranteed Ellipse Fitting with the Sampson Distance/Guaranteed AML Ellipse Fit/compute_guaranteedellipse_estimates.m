%   Function: compute_guaranteedellipse_estimates
%
%   This function is a wrapper for an approximate maximum likelihood guaranteed
%   ellipse fit due to
%
%   Z.Szpak, W. Chojnacki and A. van den Hengel
%   "Guaranteed Ellipse Fitting with the Sampson Distance"
%   Proc. 12th European Conference on Computer Vision. ECCV 
%   Firenze,Italy, oct, 2012
%
%   It first shifts all the data points to a new coordinate system so that
%   the origin of the coordinate system is at the center of the data points, 
%   and then scales all data points so that they lie more or less within
%   a unit box. It is within this transformed coordinate system that the 
%   ellipse is estimated. Since the approximate maximum likelihood method
%   is an iterative scheme, it requires an initial seed for the parameters.
%   The seed is taken to be the direct ellipse fit due to 
%
%   R. Halif and J. Flusser
%   "Numerically stable direct least squares fitting of ellipses"
%   Proc. 6th International Conference in Central Europe on Computer Graphics and Visualization. WSCG '98
%   Czech Republic,125--132, feb, 1998
%
%   Their work is a numerically stable version of
%
%   A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%   "Direct Least Squares Fitting of Ellipses"
%   IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%   The final ellipse parameters estimates are transformed 
%   back to the original data space. 
%
%   Parameters:
%
%      dataPts    - a 2xN matrix where N is the number of data points	
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
%    compute_directellipse_estimates
%    guaranteedEllipseFit
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 25/7/2012 
function [theta] = compute_guaranteedellipse_estimates(data_points)

x = data_points;
n = length(x);

% Need homogenous coordinates for normalize_transform function
x = [ x; ones( 1, n ) ];

% hartley normalise the data
[x, T] = normalise2dpts(x);

% Flusser Fit
%-----------------------------------------------------
normalised_data = x;
theta = direct_ellipse_fit(normalised_data); 
theta = theta / norm(theta);
%-----------------------------------------------------

theta = guaranteedEllipseFit( theta, x );
theta = theta/norm(theta);

a = theta(1); b = theta(2) ; c = theta(3) ; d = theta(4) ; e = theta(5); f = theta(6);
C = [a b/2 d/2 ; b/2 c e/2; d/2 e/2 f];

% denormalise C
C = T'*C*T;
aa = C(1,1);
bb = C(1,2)*2;
dd = C(1,3)*2;
cc = C(2,2);
ee = C(2,3)*2;
ff = C(3,3);
theta = [aa bb cc dd ee ff]';
theta = theta / norm(theta);


end

 



