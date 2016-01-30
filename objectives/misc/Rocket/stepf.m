%stepf.m
function step=stepf(x,s,ss)
% This is a step function stepf(x,s,ss)
% x is the input array, s the value where the step occurs
% if ss=1 the step occurs at 's' from 1 to 0 for all x
% if ss=2 the step occurs at 's' from 0 to 1 for all x
% if ss is neither 1 nor 2 the function returns the value of zero
if ss==1
step=1./(1+exp(100*(x-s)));
elseif ss==2
        step=1-1./(1+exp(100*(x-s)));
    else
        step=0.0;
end
