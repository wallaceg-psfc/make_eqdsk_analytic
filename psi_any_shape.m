function [psival] = psi_any_shape(z,y)

% Clean version
% This function returns psi normalized between 0 and 1. It will work for any of
% the equiliibrium types we have considered.
% This is meant to replace all previous functions. We may want a 2D
% option??
% psisign also needs to have been defined before this function is called.
% This works for both the pedestal and the smooth solutions. Keep in mind
% that in the pedestal case this calculates psi_hat.

global sfull cfull scfull ssfull clist sl anC bnC anS bnS nseries  powers ...
    powers_m1 hvec kvec psisign eq_option eps psimax twoD pedestal_option


nz = length(z);
ny = length(y);
%z = z*3;

if((twoD==1)&(nz>1)&(ny>1))
    % Fill a matrix
    
    z = repmat(z,ny,1);
    y = repmat(y,nz,1);
    y = y';
    z = reshape(z,ny*nz,1);
    y = reshape(y,ny*nz,1);
    nz_true = nz;
    ny_true = ny;
    nz = ny*nz;
end

psival = NaN(nz,1);

for i = 1:nz
    
    zloc = z(i);
    yloc = y(i);
    
    Xcvalue = NaN(1,4);
    %Only use one point at a time in this version,
    % but leave the array for future possible ones.
    Xsvalue = Xcvalue;
    
    zs = zloc.^powers;
    zss = repmat(zs(:),1,4);
    
    Xcvalue(:,:) = dot(zss,anC).*cos(kvec*zloc) + ...
        dot(zss,bnC).*sin(kvec*zloc);
    
    Xsvalue(:,:) = dot(zss,anS).*cos(kvec*zloc) + ...
        dot(zss,bnS).*sin(kvec*zloc);
    
    % This should be the only thing that needs to be changed between
    % the various equilibrium options
    if((eq_option == 1)|(eq_option == 2))
        % Miller profile or Double Null equilibrium
        %             res(i,j) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue);
        res(i) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue);
    elseif(eq_option == 3)
        % Single Null equilibrium
        
        %             res(i,j) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue)+...
        res(i) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue)+...
            dot(sin(hvec(1:3)*yloc),ssfull.*Xsvalue(1:3)+scfull.*Xcvalue(1:3));
        
    else
        % Something went wrong
        res = 0;
    end
    
    %         res(i,j) = res(i,j)*psisign;
    
    %     end
end

if(pedestal_option==0)
    psival = psisign*res./psimax;
elseif(pedestal_option==1)
    psival = res;
end

if((twoD==1)&(nz>1)&(ny>1))
    % Turn the 1D array into a matrix
    
    psival = reshape(psival,ny_true,nz_true);
    
end

end