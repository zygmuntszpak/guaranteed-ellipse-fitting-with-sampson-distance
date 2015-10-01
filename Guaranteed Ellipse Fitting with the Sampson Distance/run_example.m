%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the guaranteed ellipse fitting
% method described in the paper:
%
%   Z.Szpak, W. Chojnacki and A. van den Hengel
%   "Guaranteed Ellipse Fitting with the Sampson Distance"
%   Proc. 12th European Conference on Computer Vision. ECCV 
%   Firenze,Italy, oct, 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% noisy data points sampled from an ellipse
data_points = [554.4255448860487, 252.25695182609456;
    563.2198518844277, 266.393133656879;
    529.4069191555953,  293.8378723522843;
    540.4542393161609, 305.7810697735044;
    529.979807146737, 316.8126082331031;
    504.904046118169, 337.4825804687354;
    481.8704189162051, 364.2838371680541;
    444.52860559075106, 363.1246963316071;
    399.2429477714536, 377.5681669353429;
    389.3711948977463, 385.5186736607949;
    349.4010115271496,386.2628649432312;
    296.5840933162226, 399.73317608469233;
    255.11250601488572, 395.1411289783602;
    231.89205825187307, 381.61074589981473;
    182.90013818617746,384.047826837435;
    144.59336259440926,391.4128662876713;
    97.60448452233265, 391.35548699175257;
    76.28734052405164,369.8323364160187;
    25.345944215923343,367.78769391856787;
    -2.5177280792243053,339.56485858969074;
    -4.3254400982999766,305.06389938514843;
    -29.074863473090797,319.05942892395166;
    -48.73339406796903, 288.8232699112716;
    -85.21858392682377,262.78052469826156;
    -77.7559980456105,255.63059999535096;
    -63.96460001891381,240.5558431574499;
    -47.50679886929984, 221.79850196104366;
    -67.29391840704433, 209.6724162607794;
    -23.86911438928088, 176.9721535710241;
    -17.624006835223113, 174.24367926036066;
    31.38415275933525, 142.62191100176753;
    46.11460906558086, 128.39144744279602;
    65.38963556066517,126.90848784795256;
    116.43356672461559,125.50674998521659;
    135.24857840307698, 109.847540360864;
    212.24483407970993, 105.53786425855019;
    227.26351500381813, 101.51465346090531;
    255.18532643477042, 107.83975623109346;
    306.56565641797596, 100.03580557868527;
    355.07952555897765, 143.5543623673256;
    376.5664452416278, 116.9383797341007;
    416.266045280784, 102.4712728975143;
    451.06342651230636, 144.36320327683816;
    475.96962498992485,  153.3519382363760;
    507.98479900085863, 186.347765271335;
    522.3993024410767, 181.1256514766796;
    544.0355186255582, 222.30892095626606;
    561.1263519611979, 205.47120976368097;
    581.8305025046524,241.94356048025762;
    546.5064585867074, 232.51223747541997]';


%%%%%%%%%%%%%%%%%%%% Example with ALL data points %%%%%%%%%%%%%%%%%%%%%%%%
% An example of fitting to all the data points
[theta_dir]  = compute_directellipse_estimates(data_points);
[theta_guaranteed] = compute_guaranteedellipse_estimates(data_points);

% plot the data points
x = data_points';
n = length(x);
figure('Color',[1 1 1])
plot(x(:,1),x(:,2),'b.');
hold on

% determine data range
minX = min(min(x(:,1))) - 20;
minY = min(min(x(:,2))) - 20;
maxX = max(max(x(:,1))) + 20;
maxY = max(max(x(:,2)))  + 20;

% plot the direct ellipse fit
a = theta_dir(1); b = theta_dir(2); c = theta_dir(3);
d = theta_dir(4); e = theta_dir(5); f = theta_dir(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);


% plot the guaranteed ellipse fit
a = theta_guaranteed(1); b = theta_guaranteed(2); c = theta_guaranteed(3);
d = theta_guaranteed(4); e = theta_guaranteed(5); f = theta_guaranteed(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);
set(h, 'Color', [1 0 1]);
legend('DATA POINT','DIRECT ELLIPSE FIT','GUARANTEED ELLIPSE FIT','Location','SouthOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%% Example with only a portion of data points %%%%%%%%%%%%%%%%%%%%%%%%
data_points_portion = data_points(:,1:1:(end/2));
% An example of fitting to a portion of the data points
[theta_dir]  = compute_directellipse_estimates(data_points_portion);
[theta_guaranteed] = compute_guaranteedellipse_estimates(data_points_portion);

% plot the data points
x = data_points';
xportion = data_points_portion';
n = length(x);
figure('Color',[1 1 1])
%plot(x(:,1),x(:,2),'b.');
plot(xportion(:,1),xportion(:,2),'b.');
hold on

% determine data range
minX = min(min(x(:,1))) - 20;
minY = min(min(x(:,2))) - 20;
maxX = max(max(x(:,1))) + 20;
maxY = max(max(x(:,2)))  + 20;

% plot the direct ellipse fit
a = theta_dir(1); b = theta_dir(2); c = theta_dir(3);
d = theta_dir(4); e = theta_dir(5); f = theta_dir(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);


% plot the guaranteed ellipse fit
a = theta_guaranteed(1); b = theta_guaranteed(2); c = theta_guaranteed(3);
d = theta_guaranteed(4); e = theta_guaranteed(5); f = theta_guaranteed(6);
fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
h = ezplot(fh,[minX maxX minY maxY]);
axis([minX maxX minY maxY]);
set(h, 'Color', [1 0 1]);
legend('DATA POINT','DIRECT ELLIPSE FIT','GUARANTEED ELLIPSE FIT','Location','SouthOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

