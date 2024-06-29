function W_r = WindUncertainty(W_i, RESA, RESB)
%WINDUNCERTAINTY randomizes wind according to user inputs. Percent standard
%deviation is 12.5%
%
%Inputs: Initial wind vector (wind vs height), resolution of output heights vector,
%resolution of randomization (lower is not better!)
%
%Outputs: Varied wind vector


%NOTE: THIS IS DONE A LIL JANKY ATM, AREA FOR IMPROVMENT FS

%Select only wind values in your resolution
h = W_i(:, 1);

RESB = RESB/RESA;

%Vector allocation for randomization
W_u = W_i(1:RESB:end, :);

%Randomizes every RESB meters
sigma_mult = 0.15;
sigma = sigma_mult.*randn(length(W_u), 2);

%Applies randomization to new vector
W_s = zeros(length(W_i), 3);
j = 1;
for i = 1:RESB:length(W_i)
    W_s(i, 2:3) = sigma(j, :);
    j = j+1;
end

%For each point between randomizations, take weighted average of two
%randomized values to form a smooth curve
%TBH i kinda forgor how i made this work but trust
for i = 1:length(W_s)

    %Index of lower randomized value
    IL = 1+floor(i/RESB)*RESB;
    IH = 1+ceil(i/RESB)*RESB;

    %Index of higher randomized value
    WL = 1-abs(IL-i)/RESB;
    WH = 1-abs(IH-i)/RESB;
    if(IL > length(W_s)), IL = length(W_s); end
    if(IH > length(W_s)), IH = length(W_s); end

    %Value of lower and higher values
    VL = W_s(IL, :);
    VH = W_s(IH, :);

    %Weighted average (forms uncertainty curve)
    W_s2(i, :) = (VL*WL+VH*WH)/(WL+WH);

end

%Add uncertainty to wind!
W_s2 = W_s2+1;
W_r = W_i.*W_s2;

W_r(:, 2) = smoothdata(W_r(:, 2));
W_r(:, 3) = smoothdata(W_r(:, 3));


end