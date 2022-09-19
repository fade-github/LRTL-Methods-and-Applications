function [seed] = setSeed(varargin)
%SETSEED Summary of this function goes here
%   Detailed explanation goes here
params = inputParser;
params.addParameter('seed',nan,@isscalar)
params.parse(varargin{:})

seed = params.Results.seed;
if isnan(seed)
    seed = rng('shuffle');
    seed = seed.Seed;
else
    rng(seed)
end

