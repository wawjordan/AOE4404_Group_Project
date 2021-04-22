classdef inviscid_burgers1D < soln1D
    properties
        
    end
    properties (Dependent = true, SetAccess = private)
        
    end
    methods
        function this = inviscid_burgers1D(xi,inputs,varargin)
            this = this@soln1D(xi,inputs);
            defaultExactSolutionType = '';
            expectedSolutions = {};
            p = inputParser;
            valid_func = @(x) isa(x,'function_handle');
            addParameter(p,'ExactSolutionType',defaultExactSolutionType,...
               @(x) any(validatestring(x,expectedSolutions)));
            parse(p,varargin{:});
           
        end
    end
end