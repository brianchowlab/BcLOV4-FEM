function Elem = lagrange_deg2_dim1(Type_STR)
%lagrange_deg2_dim1
%
%    This defines a reference finite element to be used by FELICITY.
%    INPUT: string = 'CG'/'DG' (continuous/DIScontinuous galerkin).
%
%    Lagrange H^1 Finite Element of degree = 2, in dimension = 1.
%
%    Reference Domain: unit interval [0, 1].
%
%            |------------|--> x
%            0            1

% Copyright (c) 04-Apr-2018,  Shawn W. Walker

% name it!
Elem.Name = mfilename;

% determine continuous or discontinuous galerkin space
if (nargin==0)
    Type_STR = 'CG'; % default
end
if strcmpi(Type_STR,'cg')
    Elem.Type = 'CG';
elseif strcmpi(Type_STR,'dg')
    Elem.Type = 'DG';
else
    error('Invalid input type!  Must be ''CG'' or ''DG''.');
end

% topological dimension and domain
Elem.Dim = 1;
Elem.Domain = 'interval';

% nodal basis function definitions (in the order of their local Degree-of-Freedom index)
Elem.Basis(3).Func = []; % init
Elem.Basis(1).Func = {'2*x^2 - 3*x + 1'};
Elem.Basis(2).Func = {'x*(2*x - 1)'};
Elem.Basis(3).Func = {'-4*x*(x - 1)'};
% local mapping transformation to use
Elem.Transformation = 'H1_Trans';
Elem.Degree = 2;

% nodal variable data (in the order of their local Degree-of-Freedom index)
% (the nodal point is given in barycentric coordinates)
Elem.Nodal_Var(3).Data = []; % init
Elem.Nodal_Var(1).Data = {[1, 0], 'eval_vertex', 1};
Elem.Nodal_Var(2).Data = {[0, 1], 'eval_vertex', 2};
Elem.Nodal_Var(3).Data = {[1/2, 1/2], 'eval_cell', 1};

%%%%% topological arrangement of nodal points %%%%%
% note: multiple sets of nodal variables can be associated with each topological entity.
%       This is done by defining multiple matrix arrays within each matlab cell.
%       See the FELICITY manual for more info.

% nodes attached to vertices
Elem.Nodal_Top.V = {...
                   [1;
                    2]...
                   };
%

% nodes attached to edges
Elem.Nodal_Top.E = {...
                   [3]...
                   };
%

% nodes attached to triangles (faces)
Elem.Nodal_Top.F = {[]};

% nodes attached to tetrahedra (cells)
Elem.Nodal_Top.T = {[]};

end