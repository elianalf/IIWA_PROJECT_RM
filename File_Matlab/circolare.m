function [pe, ped, pedd] = circolare ( sc, scd, scdd, c, ro, pe, ped, pedd, l, S_new)
        
        len=length(sc)-1;
        temp = zeros(length(S_new),3);

        pc = c + [ro*cos(sc/ro)     ro*sin(sc/ro)     zeros(length(sc),1)];
        
        pcd = [-1*scd.*sin(sc/ro) ...
               scd.*cos(sc/ro) ...
               zeros(length(sc),1)];
        
        pcdd = [-1*scdd.*sin(sc/ro)-((scd).^2).*cos(sc/ro)/ro ...
                scdd.*cos(sc/ro)- ((scd).^2).*sin(sc/ro)/ro ...
                zeros(length(sc),1)];

        temp(l-len:l,:)=pc;
        temp (l+1:end,:) = temp (l+1:end,:) + pc(end,:);
        pe(l-len:end,:) = temp(l-len:end,:);

        ped(l-len:l,:)=pcd;
        pedd(l-len:l,:)=pcdd;




end