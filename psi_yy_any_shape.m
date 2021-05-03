function [psival] = psi_yy_any_shape(z,y)

% This function returns psi_yy for any of the equiliibrium types we have considered.
% This is meant to be called only for postptocessing.
% psisign also needs to have been defined before this function is called.


global sfull cfull scfull ssfull clist sl anC bnC anS bnS nseries  powers ...
    powers_m1 hvec kvec psisign eq_option eps psimax

[n1 n2] = size(z);

res = NaN(n1,n2);

for i = 1:n1
    for j = 1:n2

        zloc = z(i,j);
        yloc = y(i,j);

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
            res(i,j) = -dot(hvec.^2.*cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue);
        elseif(eq_option == 3)
            % Single Null equilibrium
    
            res(i,j) = -dot(hvec.^2.*cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue) - ...
            dot(hvec(1:3).^2.*sin(hvec(1:3)*yloc),ssfull.*Xsvalue(1:3)+scfull.*Xcvalue(1:3));

        else
            % Something went wrong
            res = 0;
        end
        
        res(i,j) = res(i,j)*psisign;
        
    end
end

psival = res./psimax;

end