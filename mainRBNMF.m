clear
warning off

q2double = @(X) double2q(X,'inverse');

GT = double(imread("1.png"));  % groundtruth

Y = double2q(GT);

[n1, n2, n3] = size(GT);
M = zeros(n1, n2);
%%%Inpaithing
rate = 0.25;    % sampling rate
omega = rand(n1 * n2 , 1) < rate;
M(omega) = 1;
omega1 = find(M);


%% algorithm for image restoration

imgN = double2q(zeros(n1,n2,n3));
imgN(omega1) = Y(omega1); 
psnr1 = psnr(double(GT)./255, q2double(imgN)./255);
ssim1 = ssim(double(GT)./255, q2double(imgN)./255);
fprintf('PSNR is %2.2f dB\n', psnr1);
fprintf('SSIM is %2.4f dB\n', ssim1);
% figure;
X1 = RBNMF(imgN, M);
% X1 = RBWNNM(imgN, M);
X1(omega1)  =imgN(omega1);
X1 = q2double(X1);
psnr1 = psnr(X1./255,double(GT)./255);
SSIM1 = ssim(double(GT)./255, X1./255);
fprintf('PSNR achieved by RBNMF is %2.2f dB\n', psnr1);
fprintf('SSIM achieved by RBNMF is %2.4f dB\n', SSIM1);

