function FS = Prune_BasisFunc(obj,FS,FM)
%Prune_BasisFunc
%
%   This removes the basis functions that are not used by the FM object.

% Copyright (c) 06-24-2012,  Shawn W. Walker

for ii = 1:length(FS.Integration)
    % find the FE spaces that are needed to evaluate contributions to the FE
    %      matrices on the current integration domain
    DoI = FS.Integration(ii).DoI_Geom.Domain.Integration_Domain;
    Coef_Funcs = FS.Integration(ii).CoefFunc.values;
    Space_Names = get_spaces_for_DoI(Coef_Funcs,FM,DoI);
    
    BF_Set = FS.Integration(ii).BasisFunc.values;
    for bi = 1:length(BF_Set)
        Current_BF = BF_Set{bi};
        
        % search through matrices to see if it is used
        USED = verify_basisfunc_is_used(Space_Names,Current_BF);
        % if it is NOT used, then remove it!
        if ~USED
            remove(FS.Integration(ii).BasisFunc, Current_BF.Space_Name);  
        end
    end
end

end

function Space_Names = get_spaces_for_DoI(Coef_Funcs,FM,DoI)

% get matrices that have a contribution from the current integration domain
FM_Int_Index = FM.Get_Integration_Index(DoI);
Matrix = FM.Integration(FM_Int_Index).Matrix.values;

Map_Names = containers.Map; % init
for mi = 1:length(Matrix)
    if ~isempty(Matrix{mi}.row_func.Space_Name)
        Map_Names(Matrix{mi}.row_func.Space_Name) = Matrix{mi}.row_func.Space_Name;
    end
    if ~isempty(Matrix{mi}.col_func.Space_Name)
        Map_Names(Matrix{mi}.col_func.Space_Name) = Matrix{mi}.col_func.Space_Name;
    end
end

% add in any spaces needed for coefficient functions
for ci = 1:length(Coef_Funcs)
    Map_Names(Coef_Funcs{ci}.Space_Name) = Coef_Funcs{ci}.Space_Name;
end

Space_Names = Map_Names.keys; % just need a cell array of names

end

function USED = verify_basisfunc_is_used(Space_Names,BasisFunc)

USED = false; % init
%disp(['Current BasisFunc Space: ', BasisFunc.Space_Name]);

for ind = 1:length(Space_Names)
    if strcmp(BasisFunc.Space_Name,Space_Names{ind})
        USED = true;
        break;
    end
end

end