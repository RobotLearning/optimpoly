% construct the IIR tap weights
function [B, A] = myButter2ndOrder(W)

    V  = tan(W * pi/2);  
    Sg = V ^ 2;
    Sp = V * [-1-1i, -1+1i] / sqrt(2);
    % Bilinear transform:
    P = (1 + Sp) ./ (1 - Sp);
    G = real(Sg / prod(1 - Sp));
    % From Zeros, Poles and Gain to numerator and denominator:
    B = G * [1, 2, 1];
    A = real(poly(P));

end