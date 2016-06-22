function OU = quantile_loss(Signal, Error, Model, t)
%QUANTILE_LOSS: The quantile loss, defined below, is used in models for predicting typical values. It can be considered as a skewed version of the absolute loss.
if nargin == 1
    OP = Signal;
    t = 3;
elseif nargin == 2;
    OP = Signal;
    t = Error;
elseif nargin == 3
        t = 5;
        OP = (Signal-Model)./Error;
end

OU = zeros(size(OP));
OU(OP>=0) = t*OP(OP>=0);
OU(OP<0)  = (1-t)*OP(OP<0);

end
