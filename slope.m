%% Slope method

function [alpha_ret] = slope(z,B)


% n=numel(B);
res=(z(2)-z(1));

alpha_ret = zeros(1,length(B));

ii = find(B <= 0.0);

alpha_ret(1:end-1)=-0.5.*(diff(log(B))./res);
alpha_ret(ii)=0.0;


end