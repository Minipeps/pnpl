function [sols,stats,fail_flag_av] = solver_pfold(C0,settings) %#codegen

cutoff_threshold = 1e12;

C0 = C0(:,settings.reorder);

C0 = bsxfun(@rdivide,C0,sqrt(sum(C0.^2,2)));


%I = settings.I;
%J = settings.J;
mon = settings.mon;
% keyboard
neq  = size(C0,1);
nmon = size(mon,2);

T = settings.T;
%T0 = settings.Tid;

C1 = C0';
vv = (C1(T));
C=(sparse(settings.II,settings.JJ,vv));

% T = settings.T;
% C = zeros(max(I{end}),nmon);
% T0 = settings.Tid;
% C1 = C0';
% C(T0) = (C1(T));
% C = sparse(C);


P = settings.P;  % permissible
R = settings.R;  % reducible
E = settings.E;  % excessive

p = neq; %# of unknonws : assume # of unknown = # of eqs
n = nmon; %# of monomials in the template
r = settings.dim; %dim of p-fold solutions
ind = settings.ind; % monomial indices action_variable * permissibles

stats.neq = size(C,1);
stats.nmon = size(C,2);
stats.basis = 0;

%
% construct modulo matrix
%
% reorder
ne = length(E);
nr = length(R);
np = length(P);

V = [E R P];
C = C(:, [E, R, P]);
% eliminate excessive monomials (done by lu in the qr paper.
% this is more general)

% the row_echelon method was used before, but the qr version is almost
% always much faster.
if(~isempty(E))
    [qq, rr, ~] = qr(full(C(:, 1 : length(E))));
    %     CC = sparse(C(:, length(E) + 1 : end));
    C2 = [rr qq'*C(:, length(E) + 1 : end)];
    kk = abs(rr(1))./abs(diag(rr)) < cutoff_threshold;
    k = find(diff(kk) == -1);
    if(isempty(k))
        k = length(kk);
    end
else
    C2 = C;
    k = 0;
end
k = k(1);
% partition C into R- and P-parts
CR = C2(k + 1 : end, ne + 1 : ne + nr);
CP = C2(k + 1 : end, end - np + 1 : end);
mm = size(CR, 1);
mx = size(CR, 1);
nn = size(CR, 2) + size(CP, 2);
if(nn - mm > r)
    error('not enough equations for that solution dimension');
end

% eliminate R-monomials (this step is included in the lu factorization
% in the paper. qr is slightly slower but more stable).
[q2, UR2_0] = qr(full(CR));
CP = q2'*CP;

% [LL,UU,PP] = lu(C(:,1:(ne+nr)));


% select basis (qr + column pivoting)
CP2 = CP(1 : nr, :);
CP3 = CP(nr + 1 : end, :);
[~, r3, e] = qr(CP3, 0);
CP4 = CP2(:, e(1 : end - r));
CB1 = CP2(:, e(end - r + 1 : end));
UP3 = r3(1 : np - r, 1 : np - r);
CB2 = r3(1 : np - r, end - r + 1 : end);
if(isempty(CP4)), CP4 = []; end
if(isempty(UP3)), UP3 = []; end
ee = [1 : ne + nr e + ne + nr];

if debug
    stats.basis = P(e(end-r+1:end));
end

V = V(ee);
% mon = mon(ee);

%
% elimination step
%
Celim = [UR2_0(1 : nr + np - r, :) [CP4; UP3]];


T = - Celim \ [CB1; CB2];
if (rcond(Celim) < 5e-18)
    fail_flag_av = 1;
else
    fail_flag_av = 0;
end

% modulo matrix
modM = zeros(r, n);
modM(:, end - r + 1 : end) = eye(r);
modM(:, ne + 1 : end - r) = T';


e = V;

%
% construct action matrix
%
% m = construct_actionmatrix(mon, modM, settings.dim, P, x);
ind2 = zeros(np,1);
ind2(P) = ind;

% [nouse,ind2] = ismember(ind2(e(n-r+1:n)),e);

ind2 = find_id(e,ind2(e(n-r+1:n)));

% ind2 = find_id(e,[ind2(e(n-r+1:n));settings.ids_a1;settings.ids_a3;settings.ids_abc]');
%
% ids_vv = ind2(end-12+1:end);
% ind2 = ind2(1:end-12);


M = zeros(n, r);
M(ind2, :) = eye(r);
m = modM * M;

[vv, ~] = eig(m');
%d=diag(dk);


%% Extract possible solutions

%setup matrix
okmon = sum(modM,1)~=0;
% mon = mon(:,e);
% MM = mon(:,okmon)';

%setup vv
vv = modM(:,okmon)'*vv;

% PRid = e(end-nr-np+1:end);
% PRid = PRid(okmon-ne);


% ids_abc = find_id(e,settings.ids_abc)'- ne;
mmid = settings.p_inverter;

%     ids_ac2_list = find_id(e,settings.ids_ac2) - ne;

ids_a13_list = find_id(e,[settings.ids_a1;settings.ids_a3]) - ne;

%     ids_a2c_list = find_id(e,[settings.ids_a2c]) - ne;

ids_a1       = ids_a13_list(1:4);

sols = [];
sid_rl = cell(1:p);
idsl = cell(1:p);
vv1l = cell(1:p);
vv2l = cell(1:p);
vv1rl = cell(1:p);
constant = cell(1:p);
mm = zeros(1:p);
mx = zeros(1:p);
for ii = 1:p
    ids_a13 = ids_a13_list([ii,p+ii]);

    sid  = ii;
    sid_rl{ii} = 1:p;
    sid_rl{ii}(sid)=[];

    % 
    ids  = [ids_a13];

    idsl{ii} = ids;

    % a3/a1
    vv1l{ii} = sqrt(vv(ids(2),:)./vv(ids(1),:));


    realid = imag(vv1l{ii}) == 0;
    vv1rl{ii} = vv1l{ii}(:,realid);


    constant{ii} = vv(ids(1),:)./vv1l{ii};


    fail_flag = isempty(vv1rl{ii});
    if ~fail_flag && min(abs(real((vv1l{ii})))) > 1e-7
        mm(ii) = min(norm(vv1l{ii}));
        mx(ii) = max(abs(vv1rl{ii}));
    else
        mm(ii) = inf;
        mx(ii) = 1;
    end


    % %             ids_a2c = ids_a2c_list(3*(ii-1)+(1:3));
    vv2l{ii}  =(vv(ids_a1(sid_rl{ii}),realid)./(ones(3,1)*(constant{ii}(:,realid))));
    solsl = [ [vv1rl{ii} -vv1rl{ii}]; [vv2l{ii} -vv2l{ii}] ] ;
    solsl([sid,sid_rl{ii}],:) = solsl;
    sols = [sols solsl];
end

realid = (sum(imag(sols)<1e-4)==4);
sols   = sols(:,realid);

end % function end

function ids_reorder = find_id (list,ids)
% slightly faster than ismember
ccc = zeros(1,length(list));
ccc(ids) = 1:length(ids);
ccc      = ccc(list);
[nouse,ids_reorder,ff] = find(ccc);
ids_reorder(ff) = ids_reorder;

end
