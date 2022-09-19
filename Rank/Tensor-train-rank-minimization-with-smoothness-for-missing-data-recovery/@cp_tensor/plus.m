function C = plus(A,B)
%PLUS Binary addition for cp_tensor.
%
%   C = plus(A,B) adds to CP tensors of the same size, and the
%   result is a CP_tensor of the same size.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.


if (isa(A,'cp_tensor') && isa(B,'cp_tensor'))    

    if ~( issamesize(A,B) )
	error('Tensor size mismatch.')
    end

    lambda = [A.lambda; B.lambda];    
    M = order(A);
    for m = 1 : M
        u{m} = [A.u{m} B.u{m}];
    end 
    C = cp_tensor(lambda, u);
    return;
end

error('Use plus(full(A),full(B)).');
