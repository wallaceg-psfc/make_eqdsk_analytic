function [psival] = psi_xx_any_shape(x,y)

% This function returns psi_zz for any of the equiliibrium types we have considered.
% This is meant to be called only for postptocessing.
% psisign also needs to have been defined before this function is called.

global sfull cfull scfull ssfull clist sl anC bnC anS bnS nseries  powers ...
    powers_m1 hvec kvec psisign eq_option eps psimax alphasol nu epsbar 

[n1 n2] = size(x);

res = NaN(n1,n2);

lamtest = alphasol*sqrt(epsbar*nu);

for i = 1:n1
    for j = 1:n2

        xloc = x(i,j);
        yloc = y(i,j);

        Xcvalue = NaN(1,4);
        %Only use one point at a time in this version,
        % but leave the array for future possible ones.
        Xsvalue = Xcvalue;

        xs = xloc.^powers;
        xss = repmat(xs(:),1,4);

        Xcvalue(:,:) = dot(xss,anC).*cos(kvec*xloc) + ...
        dot(xss,bnC).*sin(kvec*xloc);

        Xsvalue(:,:) = dot(xss,anS).*cos(kvec*xloc) + ...
        dot(xss,bnS).*sin(kvec*xloc);
    
        % Since we need a second derivative, we just redefine the XC and XS
        % values to their second derivatives.

        d2_term = -(kvec.^2+lamtest^2*xloc)/(1+epsbar*xloc);
        Xcvalue = Xcvalue .* d2_term;
        Xsvalue = Xsvalue .* d2_term;
        
        % This should be the only thing that needs to be changed between
        % the various equilibrium options
        if((eq_option == 1)|(eq_option == 2))
            % Miller profile or Double Null equilibrium
            res(i,j) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue);
        elseif(eq_option == 3)
            % Single Null equilibrium
    
            res(i,j) = dot(cos(hvec*yloc),sfull.*Xsvalue+cfull.*Xcvalue)+...
            dot(sin(hvec(1:3)*yloc),ssfull.*Xsvalue(1:3)+scfull.*Xcvalue(1:3));

        else
            % Something went wrong
            res = 0;
        end
        
        res(i,j) = res(i,j)*psisign;
        
    end
end

psival = res./psimax;

end