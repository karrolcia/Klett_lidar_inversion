%% Klett inversion procedure


function [alpha_ret, fac, fac2] = simple_inversion(i_o, alpha, z, B, cor1, cor2)

    
    res=(z(2)-z(1));

    n=numel(z);
    int_z=zeros(1,n);
    int_z(1) = B(1).* res;
    
%     i_o=-99;
%     for i=1:n 
%        if (i_o == -99) 
%           if (z(i) < z_o) 
%              i_o=i;
%           end
%        end
%     end

    for i=2:n 
       int_z(i)=int_z(i-1)+B(i-1).*res;
    end

    fac=zeros(1,n);
    fac2=ones(1,n);
    fac(1) = 0.5;
    fac2(1) = 1.0;
  for i=1:n
        
     if (cor1 == 0) 
          fac(:) = 0.5;
     end     
       if cor1 == 1
          if alpha(i) > 0.0 
             fac(i)=(exp(res*alpha(i))-1.0)/(exp(res*alpha(i))...
                 -exp(-1.0*res*alpha(i)));
             fac2(i)=(2*res*alpha(i))/...
                 (exp(res*alpha(i))-exp(-1.0*res*alpha(i)));
          else
             fac(i)=0.5;
             fac2(i)=1.0;
          end
      end
       
       int_z(i)=int_z(i)+B(i)*res*fac(i);
  end
 
  int_zo=int_z(i_o);
    
 int =-1.0*(int_z-int_zo);
  
  if cor2 == 0
       fac2(:)=1.0;
  end
  
  alpha_ret=B.*fac2./(fac2(i_o)*B(i_o)./alpha(i_o)+2.0.*int);

 end
       