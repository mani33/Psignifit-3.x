function ci = getCI ( inference, cut, p,varargin )
% ci = getCI ( inference, cut )
%
% Get the confidence interval from the inference object
%
% inference should be a struct as generated either by BayesInference or BootstrapInference
%
% cut should be the index of the desired cut
%
% p should be the desired coverage of the confidence interval (p<1)
%
%
% This file is part of psignifit3 for matlab (c) by Ingo Fründ
% Modified by MS 2012-10-09
args.param = 'threshold'; % 'slope' or 'threshold' for which we need the CI
args = parseVarArgs(args,varargin{:});

notin = 1-p;
probs = [0.5*notin,1-0.5*notin];

if strcmp ( inference.call, 'bootstrap' )
    bias = inference.bias_thres(cut);
    acc  = inference.acc_thres(cut);
    probs = normcdf( bias + ( norminv(probs) + bias ) ./ (1-acc*(norminv(probs) + bias )) );
end;

switch args.param
    case 'threshold'
        ci = prctile ( inference.mcthres(:,cut), 100*probs );
    case 'slope'
        ci = prctile ( inference.mcslopes(:,cut), 100*probs );
end
% ci = prctile ( inference.mcthres(:,cut), 100*probs );
