function W_r = WindUncertainty(W_i, sigma_s, sigma_o)
%WINDUNCERTAINTY randomizes wind according to user inputs. Percent standard
%deviation is 12.5%
%
%Inputs: Initial wind vector (wind vs height), resolution of output heights vector,
%resolution of randomization (lower is not better!)
%
%Outputs: Varied wind vector


%NOTE: THIS IS DONE A LIL JANKY ATM, AREA FOR IMPROVMENT FS

%Split Data
h = W_i(:, 1);
W = W_i(:, 2:3);

err_s = 1-sigma_s/100*randn(1, 1);
err_o = randn(1, 1)*sigma_o;

W = W*err_s+err_o;

gust = 15;
fs = 1/(h(2)-h(1));
fd = 1/1000;
noise = randn(length(h), 1);
noise = lowpass(noise, fd, fs)*gust;
W=W+noise;

plot(h, noise)



W_r = [h, W];


end