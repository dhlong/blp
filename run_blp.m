% % for OTC data
% A = importdata('OTC.csv');
% data = A.data;
% 
% A = importdata('OTC_RCDraws_forMatlab.csv');
% nu = A.data(:,3:end);
% 
% A = importdata('OTC_Demographics_forMatlab.csv');
% demogr = A.data(:,3:end);
% 
% 
% cdid    = data(:,3);
% share   = data(:,5);
% outshr  = data(:,6);
% x1      = data(:,7:19);
% IV      = data(:,20:end);
% 
% x2      = [ones(size(x1,1),1) x1(:,[1 3])];
% Z       = [x1(:,2:end) IV];

load('BLP_data.mat');

x1 = [const hpwt air mpg space log(price)];
x2 = [const hpwt air mpg space];

n.k1 = size(x1,2);
n.k2 = size(x2,2);
n.obs = size(x1,1);
n.cdid = numel(unique(cdid));

% instruments
firmcdid    = firmid*max(cdid) + cdid;
charid      = kron((1:n.k1)', ones(n.obs, 1));

sumfirm     = accumarray([repmat(firmcdid, n.k1, 1) charid], x1(:));
sumcdid     = accumarray([repmat(cdid, n.k1, 1) charid], x1(:));

ivfirm      = sumfirm(firmcdid,:) - x1;
ivcdid      = sumcdid(cdid,:) - x1;

Z           = [x1(:, 1:end-1) ivfirm(:, 1:end-1) ivcdid(:, 1:end-1)];

% random draws
n.draws     = 100;
n.demogr    = 0;
nu          = normrnd(0,1,[n.cdid n.draws*n.k2]);
demogr      = zeros(n.cdid, 0);

rng('default');

global countf;
countf = 0;

% x1 = [const hpwt log(price)];

% iv1 = zeros(size(x1) - [0 2]);
% for market = unique(cdid)'
%     for firm = unique(firmid)'
%         filter = market == cdid & firm == firmid;
%         iv1(filter,:) = repmat(sum(x1(filter,2:end-1),1), sum(filter), 1);
%     end
% end
% 
% iv2 = zeros(size(x1) - [0 2]);
% for market = unique(cdid)'
%     filter = market == cdid;
%     iv2(filter,:) = repmat(sum(x1(filter,2:end-1),1), sum(filter), 1);
% end
% 
% iv2 = iv2 - iv1;
% iv1 = iv1 - x1(:,2:end-1);
% Z = [x1(:,1:end-1) iv1 iv2];

ZZ = Z'*Z;
XZ = x1'*Z;
XZ_ZZ = XZ/ZZ;
Gamma = (XZ_ZZ*XZ')\XZ_ZZ;

% nu = normrnd(0, 1, [size(x2,2), 100]);
nu = reshape(nu, n.draws*n.cdid, n.k2);
demogr = reshape(demogr, n.draws*n.cdid, n.demogr);

cdidrep     = repmat(cdid, n.draws, 1);
drawidrep   = kron((1:n.draws)',ones(n.obs,1));
cdidmap     = cdidrep + (drawidrep-1)*n.cdid;
x2rep       = repmat(x2, n.draws, 1);

SigmaConstraint = eye(size(x2,2))>0;
% SigmaConstraint(3,3) = 1;
PiConstraint = [];
n.Sigma = sum(SigmaConstraint(:));
n.Pi = sum(PiConstraint(:));


dvcdid = bsxfun(@eq, sparse(sort(unique(cdid))), cdid');

data = v2struct(x1, x2, Z, ZZ, XZ, Gamma, share, outshr, SigmaConstraint, PiConstraint, ...
    cdid, cdidrep, drawidrep, cdidmap, dvcdid, nu, demogr, cdidmap, n);


Sigma = 0*ones([n.Sigma 1]);
Pi = 0*ones([n.Pi 1]);

if exist('delta.mat','file')
    delete delta.mat;
end


options = optimset(...
    'Display','iter', ...
    'MaxFunEvals', 10000);
theta = fminunc(@(theta) gmmobj(theta, data, n), [Sigma;Pi], options);


