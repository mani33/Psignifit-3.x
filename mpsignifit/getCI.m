function ci = getCI ( inference, cut, p, param, bca )
% ci = getCI ( inference, cut, p, param, bca )
%
% Get the confidence interval from the inference object
%
% inference should be a struct as generated either by BayesInference or BootstrapInference
%
% cut should be the index of the desired cut
%
% p should be the desired coverage of the confidence interval (p<1)
%
% 'param' should be 'slope' or 'threshold' for which we need the CI.
% 'bca' set to true to use Bias Corrected and Accelerated CI estimation.

% This file is part of psignifit3 for matlab (c) by Ingo Fründ
% 
% example:   getCI(fitData, 1, 0.95, 'slope')
%            getCI(fitData, 1, 0.95, 'threshold', true)
%
%
% Modified by MS 2013-01-30

% Default values
if nargin < 4
    param = 'threshold';
    bca = true;
end

if nargin < 5
    bca = true;
end


notin = 1-p;
probs = [0.5*notin,1-0.5*notin];

if bca && strcmp ( inference.call, 'bootstrap' )
    bias = inference.bias_thres(cut);
    acc  = inference.acc_thres(cut);
    probs = normcdf( bias + ( norminv(probs) + bias ) ./ (1-acc*(norminv(probs) + bias )) );
end;

switch param
    case 'threshold'
        ci = prctile ( inference.mcthres(:,cut), 100*probs );
    case 'slope'
        ci = prctile ( inference.mcslopes(:,cut), 100*probs );
    otherwise
        error('param should be either threshold or slope')
end
