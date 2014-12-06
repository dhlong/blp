function G = gmmobj(theta, data, n)
%GMMOBJ Summary of this function goes here
%   Detailed explanation goes here

global countf;

% Contraction mapping

Sigma = zeros(size(data.SigmaConstraint));
Sigma(data.SigmaConstraint>0) = theta(1:n.Sigma);

Pi = zeros(size(data.PiConstraint));
Pi(data.PiConstraint>0) = theta(n.Sigma+1:n.Sigma+n.Pi);

delta = log(data.share) - log(data.outshr);

% if exist('delta.mat','file')
%     load delta.mat;
% end

% beta_i = data.nu*Sigma' + data.demogr*Pi';
beta_i = data.nu*Sigma';
beta_i = reshape(beta_i, n.cdid, n.draws, n.k2);
beta_i = beta_i(data.cdid, :, :);

x2 = reshape(data.x2, n.obs, 1, n.k2);
mu = sum(bsxfun(@times, x2, beta_i),3);


toler = 1e-8;
converged = false;
count = 0;

while ~converged
    delta_0 = delta;
    
    
    %     mu = zeros(n.obs, n.draws);
    %     for i = 1:n.draws
    %         mu(:,i) = sum(data.x2.*beta_i(:,i:n.draws:(n.k2*n.draws)),2);
    %     end
    %
    %     mu = data.x2*Sigma*data.nu;
    s = exp(bsxfun(@plus, delta, mu));
    
    %     idx1 = repmat(data.cdid, 1, size(meanval,2));
    %     idx2 = repmat(1:size(meanval,2), size(data.cdid,1), 1);
    %     labels = [idx1(:) idx2(:)];
    %     maxmeanval = accumarray(labels,meanval(:),[],@max);
    %     normalizedmeanval = meanval - maxmeanval(data.cdid,:);
    
    ss = accumarray([data.cdidrep data.drawidrep], s(:)) + 1;
    s = s./ss(data.cdid,:);
    delta = delta + log(data.share) - log(mean(s,2));
    
    gap = max(abs(delta(:) - delta_0(:)));
    converged = gap < toler;
    count = count + 1;
    %     fprintf('   contraction mapping, loop = %d, gap = %.8f\n',count,gap);
end

save delta.mat delta
%
% c = zeros(size(data.share));
% for cdid = unique(data.cdid)
%     for firm = unique(data.firmid(filter))
%         filter = (cdid == data.cdid) & (firm == data.firmid);
%         ss = s(filter,:);
%         vv = data.nu(filter,:);
%         Delta_price = -(ss.*(alpha + vv*Sigma(:,1)))'*ss + diag(sum(ss.*(alpha + vv*Sigma(:,1)),2));
%         c(filter) = p - Delta_price\data.share(filter);
%     end
% end
%
%
% beta = data.XOXiXO*[delta;ln(c)];
%
% beta_d = beta(1:n.beta_d);
% beta_s = beta(1:n.beta_s);
%
% ksi = delta - data.x1*beta_d - alpha*data.price;
% omega = ln(c) - data.x2*beta_s;
%
% G = [ksi' omega']*data.Omega*[ksi;omega];


% demand first
beta = data.Gamma*(data.Z'*delta);
ksi = delta - data.x1*beta;
ksi_Z = ksi'*data.Z;
G = ksi_Z/data.ZZ*ksi_Z';

if countf == 50
    display(beta);
    display(diag(Sigma));
    countf = 0;
else
    countf = countf + 1;
end
end

