function OU = huber_loss(Signal, Error, Model, c)
%HUBER_LOSS: The Huber loss, defined below, is used mostly in real-valued regression, which is a smoothed version of the absolute loss.
%   See https://en.wikipedia.org/wiki/Huber_loss
if nargin == 1
    OP = Signal;
    c = 3;
elseif nargin == 2;
    OP = Signal;
    c = Error;
elseif nargin == 3
        c = 3;
        OP = (Signal-Model)./Error;
end

OU = zeros(size(OP));
OU(abs(OP)<=c) = 0.5*OP(abs(OP)<=c).^2;
OU(abs(OP)>c) = c*(abs(OP(abs(OP)>c)) - 0.5*c);

end

% ((abs(t) < c) * 0.5 * t ** 2 + (abs(t) >= c) * -c * (0.5 * c - abs(t)))