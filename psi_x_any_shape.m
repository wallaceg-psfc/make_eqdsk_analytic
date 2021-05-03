function [psival] = psi_x_any_shape(x,y)

% This function returns psi_x for any of the equiliibrium types we have considered.
% This is meant to be called only for postptocessing.
% psisign also needs to have been defined before this function is called.

global sfull cfull scfull ssfull clist sl anC bnC anS bnS nseries  powers ...
    powers_m1 hvec kvec psisign eq_option eps psimax

[n1 n2] = size(x);

res = NaN(n1,n2);

for i = 1:n1
    for j = 1:n2

        xloc = x(i,j);
        yloc = y(i,j);

        Xcvalue = NaN(1,4);
        %Only use one point at a time in this version,
        % but leave the array for future possible ones.
        Xsvalue = Xcvalue;
        Xcprimvalue = Xcvalue;
        Xsprimvalue = Xcvalue;

        xs = xloc.^powers;
        xss = repmat(xs(:),1,4);
        
        xs_m1 = powers.*xloc.^powers_m1;
        xss_m1 = repmat(xs_m1(:),1,4);

        Xcvalue(:,:) = dot(xss,anC).*cos(kvec*xloc) + ...
        dot(xss,bnC).*sin(kvec*xloc);

        Xsvalue(:,:) = dot(xss,anS).*cos(kvec*xloc) + ...
        dot(xss,bnS).*sin(kvec*xloc);

        Xcprimvalue(:,:) = dot(xss_m1,anC).*cos(kvec*xloc) + ...
            dot(xss_m1,bnC).*sin(kvec*xloc) + ...
            kvec.*(dot(xss,bnC).*cos(kvec*xloc) - ...
            dot(xss,anC).*sin(kvec*xloc));

        Xsprimvalue(:,:) = dot(xss_m1,anS).*cos(kvec*xloc) + ...
            dot(xss_m1,bnS).*sin(kvec*xloc) + ...
            kvec.*(dot(xss,bnS).*cos(kvec*xloc) - ...
            dot(xss,anS).*sin(kvec*xloc));

        % This should be the only thing that needs to be changed between
        % the various equilibrium options
        if((eq_option == 1)|(eq_option == 2))
            % Miller profile or Double Null equilibrium
            res(i,j) = dot(cos(hvec*yloc),sfull.*Xsprimvalue+cfull.*Xcprimvalue);
        elseif(eq_option == 3)
            % Single Null equilibrium
    
            res(i,j) = dot(cos(hvec*yloc),sfull.*Xsprimvalue+cfull.*Xcprimvalue)+...
            dot(sin(hvec(1:3)*yloc),ssfull.*Xsprimvalue(1:3)+scfull.*Xcprimvalue(1:3));

        else
            % Something went wrong
            res = 0;
        end
        
        res(i,j) = res(i,j)*psisign;
        
    end
end

psival = res./psimax;

end