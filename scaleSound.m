function y = scaleSound(x)
% Scales a matrix appropriately to between +/- 1 for
% wavwrite() for example.

clipbuffer = 0.95;

y = clipbuffer*x./max(abs(x(:)));