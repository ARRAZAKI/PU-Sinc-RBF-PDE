% function main_rbf_pum_rbf_switch
clc; close all;

% =========================================================
% Script RBF-PU pour : Poisson, Helmholtz, Burgers nonlineaire, probleme singulier
% Choix du noyau RBF : 'sinc' ou 'mq' (multiquadric)ou 'gaussian'
% Choix du domaine : 'square' ou 'cassini'
% =========================================================

% ===== Parametres utilisateur =====
opt.problem = 'poisson';   % 'poisson', 'helmholtz', 'burgers_nl', 'singular_linear'
opt.domain  = 'square';      % 'square' ou 'cassini'
opt.epsmode = 'loocv';        % garder loocv pour la comparaison principale
opt.epsmode_list = {'loocv'};%,'loocv_ga','loocv_cond'
opt.cond_limit = 1e12;         % seuil numerique a documenter
opt.nRepeats = 1;             
opt.data_perturbation = false; % pour les EDP lineaires uniquement

% Si domain='square'  : Nlist represente le nombre cible de noeuds.
% Si domain='cassini' : Nlist represente le nombre de subdivisions de la grille.
opt.Nlist   = [169; 255; 441; 625];   %;  nombre total cible de noeuds pour domain='square'
opt.Nlist_cassini = [7; 13; 21; 27];

% Choix du noyau RBF :
%   {'sinc'}       : Sinc-RBF seulement
%   {'mq'}         : Multiquadric seulement
%   {'sinc','mq'}  : comparaison des deux noyaux
opt.rbf_list = {'sinc'}; % 'mq','gaussian' 

% Pour Helmholtz seulement : tester plusieurs valeurs de lambda
% Pour Poisson et Burgers, le script utilise automatiquement lambda = 1
opt.lambda_list = [0.5; 1; 10; 50; 100];

opt.epsfix  = 7;
opt.epscand = linspace(1,8,10);   
% LOOCV de type c_opt pour les EDP lineaires : c = pi*epsilon
opt.c_bounds = [0.01,10];
opt.loocv_tol = 1e-3;
opt.s_d     = 60;
opt.show    = false;   % eviter les nombreuses figures pendant les repeats

% Parametres du domaine Cassini
opt.cassini_nn      = 140;
opt.cassini_xshift  = 0;      % garder Cassini original sans decalage
opt.singular_x_tol  = 1e-12;  % eviter x=0 comme point de collocation

% ===== Gestion de lambda =====
if strcmpi(opt.problem,'helmholtz')
    lambda_list = opt.lambda_list(:)';   % plusieurs lambdas possibles
else
    lambda_list = 1;                     % lambda n'est pas utilise pour poisson/burgers
end

if strcmpi(opt.domain,'cassini')
    Nlist = opt.Nlist_cassini(:)';
else
    Nlist = opt.Nlist(:)';
end

AllResults = table();

% ===== Experiences repetitions x criteres x noyaux =====
for irep = 1:opt.nRepeats
for imode = 1:numel(opt.epsmode_list)
opt.epsmode = opt.epsmode_list{imode};
% ===== Boucle sur le noyau RBF =====
for irbf = 1:numel(opt.rbf_list)

    rbf_type = lower(opt.rbf_list{irbf});
    if ~strcmp(rbf_type,'sinc') && ~strcmp(opt.epsmode,'loocv')
        continue;
    end
    check_rbf_type(rbf_type);

    fprintf('\n########################################\n');
    fprintf('Noyau RBF : %s\n', upper(rbf_type));
    fprintf('########################################\n');

    % ===== Boucle sur lambda =====
    for il = 1:numel(lambda_list)

        lambda = lambda_list(il);
        prob = get_problem(opt.problem,lambda);

        RMSE = zeros(numel(Nlist),1);
        LINF = zeros(numel(Nlist),1);
        TIME = zeros(numel(Nlist),1);
        Nactual = zeros(numel(Nlist),1);
        Ulast = [];
        EPSMEAN = nan(numel(Nlist),1);
        EPSMIN = nan(numel(Nlist),1);
        EPSMAX = nan(numel(Nlist),1);
        CONDMEAN = nan(numel(Nlist),1);
        CONDMAX = nan(numel(Nlist),1);
        BADFRAC = nan(numel(Nlist),1);
        DATAAMPL = nan(numel(Nlist),1);

        fprintf('\n==============================\n');
        fprintf('Probleme : %s | RBF = %s | lambda = %g\n', opt.problem, rbf_type, lambda);
        fprintf('==============================\n');

        % ===== Boucle sur les tailles N =====
        for k = 1:numel(Nlist)
            N = Nlist(k);

            switch lower(opt.domain)
                case 'square'
                    [data,indI,indB] = generate_nodes(N,irep);
                    bdpt = [];
                    xeval = make_grid(opt.s_d);
                case 'cassini'
                    if strcmpi(opt.problem,'singular_linear')
                        xzero_tol = opt.singular_x_tol;
                    else
                        xzero_tol = 0;
                    end
                    [data,indI,indB,bdpt] = generate_nodes_cassini( ...
                        N,opt.cassini_nn,opt.cassini_xshift,xzero_tol,irep);
                     xeval = make_grid_cassini(opt.s_d,bdpt);
                otherwise
                    error('Domaine inconnu. Utiliser ''square'' ou ''cassini''.');
            end

            Nactual(k) = size(data,1);
            m_d = max(2,floor(sqrt(Nactual(k))/4));

            tic

            PU = pu_setup(data,m_d);
            PU = prepare_local_pu(data,PU);

            % Choix du parametre epsilon patch par patch
            if ismember(lower(opt.epsmode),{'loocv','econm_c','loocv_grid','loocv_cond','loocv_ga'})
                eps_patch = zeros(PU.M,1);
                for jp = 1:PU.M
                    idx = PU.patches{jp};
                    if strcmpi(opt.problem,'burgers_nl')
                        % Premier LOOCV Burgers : linearisation initiale ubar=0,
                        % donc operateur -Delta.
                        eps_patch(jp) = choose_eps_burgers_patch( ...
                            data(idx,:),ismember(idx,indB),prob, ...
                            opt.c_bounds,opt.loocv_tol,rbf_type,opt.epsmode,opt.cond_limit);
                    else
                        eps_patch(jp) = choose_eps_pde_patch( ...
                            data(idx,:),ismember(idx,indB),prob,opt.problem,lambda, ...
                            opt.c_bounds,opt.loocv_tol,rbf_type,opt.epsmode,opt.cond_limit);
                    end
                end
            else
                eps_patch = opt.epsfix*ones(PU.M,1);
            end

            % Matrices locales RBF
            PU = attach_local_mats(PU,eps_patch,rbf_type);

            % Assemblage des operateurs PU
               if strcmpi(opt.problem,'burgers_nl') || strcmpi(opt.problem,'singular_linear')
                    [Ax,Ay,Axx,Ayy] = build_pu_operators(data,PU,eps_patch,rbf_type);
                  else
                    [Axx,Ayy] = build_pu_laplacian(data,PU,eps_patch,rbf_type);
               end

            if strcmpi(opt.problem,'singular_linear')

                xI = data(indI,1);
                if any(abs(xI) < 1e-12)
                    error('Le probleme singulier contient 1/x et 1/x^2 : eviter les points interieurs avec x=0.');
                end
            
                D1 = spdiags(1./xI,0,numel(indI),numel(indI));
                D2 = spdiags(1./(xI.^2),0,numel(indI),numel(indI));
            
                LI = Axx(indI,indI) + D2*Ayy(indI,indI) + D1*Ax(indI,indI);
                LB = Axx(indI,indB) + D2*Ayy(indI,indB) + D1*Ax(indI,indB);
              else
                LI = Axx(indI,indI) + Ayy(indI,indI);
                LB = Axx(indI,indB) + Ayy(indI,indB);
            
            end

            fI = prob.f(data(indI,1),data(indI,2));
            gB = prob.g(data(indB,1),data(indB,2));

            % ===== Resolution =====
            if strcmpi(opt.problem,'poisson')

                uI = (-LI) \ (fI + LB*gB);

            elseif strcmpi(opt.problem,'helmholtz')

                uI = (-LI + lambda*speye(numel(indI))) \ (fI + LB*gB);

            elseif strcmpi(opt.problem,'singular_linear')

                uI = LI \ (fI - LB*gB);

            elseif strcmpi(opt.problem,'burgers_nl')

                LxI = Ax(indI,indI);  LxB = Ax(indI,indB);
                LyI = Ay(indI,indI);  LyB = Ay(indI,indB);

                uI = zeros(numel(indI),1);
                maxit = 30;
                tol = 1e-10;

                % ===== Premiere resolution de Burgers =====
                % epsilon vient du premier LOOCV avec ubar=0, donc -Delta.
                for it = 1:maxit
                    ux = LxI*uI + LxB*gB;
                    uy = LyI*uI + LyB*gB;

                    F = -LI*uI + uI.*ux + uI.*uy - fI - LB*gB;

                    if norm(F,inf) < tol
                        break;
                    end

                    J = -LI + diag(ux + uy) + diag(uI)*LxI + diag(uI)*LyI;
                    du = -J\F;
                    uI = uI + du;

                    if norm(du,inf) < tol
                        break;
                    end
                end

                % ===== Une seule mise a jour LOOCV pour Burgers =====
                % Utilise uniquement la solution numerique obtenue ci-dessus.
                if ismember(lower(opt.epsmode),{'loocv','econm_c','loocv_grid','loocv_cond','loocv_ga'})
                    ubar = zeros(size(data,1),1);
                    ubar(indI) = uI;
                    ubar(indB) = gB;
                    uxbar = Ax*ubar;
                    uybar = Ay*ubar;

                    for jp = 1:PU.M
                        idx = PU.patches{jp};
                        eps_patch(jp) = choose_eps_burgers_patch( ...
                            data(idx,:),ismember(idx,indB),prob, ...
                            opt.c_bounds,opt.loocv_tol,rbf_type,opt.epsmode,opt.cond_limit, ...
                            ubar(idx),uxbar(idx),uybar(idx));
                    end

                    % Reconstruire uniquement ce qui depend de epsilon.
                    PU = attach_local_mats(PU,eps_patch,rbf_type);
                    [Ax,Ay,Axx,Ayy] = build_pu_operators(data,PU,eps_patch,rbf_type);

                    LI = Axx(indI,indI) + Ayy(indI,indI);
                    LB = Axx(indI,indB) + Ayy(indI,indB);
                    LxI = Ax(indI,indI);  LxB = Ax(indI,indB);
                    LyI = Ay(indI,indI);  LyB = Ay(indI,indB);

                    % Relancer Newton en prenant la premiere solution comme initialisation.
                    for it = 1:maxit
                        ux = LxI*uI + LxB*gB;
                        uy = LyI*uI + LyB*gB;

                        F = -LI*uI + uI.*ux + uI.*uy - fI - LB*gB;

                        if norm(F,inf) < tol
                            break;
                        end

                        J = -LI + diag(ux + uy) + diag(uI)*LxI + diag(uI)*LyI;
                        du = -J\F;
                        uI = uI + du;

                        if norm(du,inf) < tol
                            break;
                        end
                    end
                end

            else
                error('Probleme inconnu.');
            end

            u = zeros(size(data,1),1);
            u(indI) = uI;
            u(indB) = gB;

            Uh = eval_pu_solution(xeval,data,u,PU,eps_patch,rbf_type);
            TIME(k) = toc;

            err = Uh - prob.u(xeval(:,1),xeval(:,2)); % solution exacte pour EVALUATION seulement
            RMSE(k) = sqrt(mean(err.^2));
            LINF(k) = max(abs(err));
            Ulast = Uh;
            EPSMEAN(k)=mean(eps_patch);
            EPSMIN(k)=min(eps_patch);
            EPSMAX(k)=max(eps_patch);
            localrcond = nan(PU.M,1);
            for jdiag=1:PU.M
                Kraw = rbf_kernel_matrix(distmat(PU.Xloc{jdiag},PU.Xloc{jdiag}), ...
                    eps_patch(jdiag),rbf_type);
                localrcond(jdiag)=rcond(Kraw);
            end
            cv=1./localrcond;
            CONDMEAN(k)=mean(cv);
            CONDMAX(k)=max(cv);
            BADFRAC(k)=mean(cv>opt.cond_limit);
            % Sensibilite a une petite perturbation deterministe du second membre.
            % Ne modifie ni les epsilon ni la solution officielle Uh.
            if opt.data_perturbation && ~strcmpi(opt.problem,'burgers_nl')
                if strcmpi(opt.problem,'poisson')
                    Ag = -LI; bg = fI + LB*gB;
                elseif strcmpi(opt.problem,'helmholtz')
                    Ag = -LI+lambda*speye(numel(indI)); bg = fI+LB*gB;
                else
                    Ag = LI; bg = fI-LB*gB;
                end
                perturb=1e-8*max(norm(bg,inf),1)*sin((1:numel(bg))');
                du_test=Ag\(bg+perturb)-uI;
                relSol=norm(du_test,2)/max(norm(uI,2),eps);
                relRhs=norm(perturb,2)/max(norm(bg,2),eps);
                DATAAMPL(k)=relSol/max(relRhs,eps);
            end

            fprintf('RBF=%s | lambda=%g | Ninput=%d | Nactual=%d | RMSE=%.3e | LINF=%.3e | TIME=%.2fs\n', ...
                rbf_type, lambda, N, Nactual(k), RMSE(k), LINF(k), TIME(k));
        end

        RBF = repmat({rbf_type},numel(Nlist),1);
        T = table(RBF,repmat({opt.domain},numel(Nlist),1), ...
            repmat({opt.epsmode},numel(Nlist),1),irep*ones(numel(Nlist),1), ...
            lambda*ones(numel(Nlist),1),Nlist(:),Nactual,RMSE,LINF,TIME, ...
            EPSMEAN,EPSMIN,EPSMAX,CONDMEAN,CONDMAX,BADFRAC,DATAAMPL, ...
            'VariableNames',{'RBF','Domain','Mode','Repeat','lambda','Ninput','N', ...
            'RMSE','LINF','Time_s','Eps_mean','Eps_min','Eps_max', ...
            'Cond_mean','Cond_max', ...
            'Bad_fraction','Data_amplification'});

        fprintf('\n=== Resultats pour RBF = %s | lambda = %g ===\n', rbf_type, lambda);
        disp(T);

        AllResults = [AllResults; T]; 

        if opt.show
            figure; loglog(Nactual,RMSE,'-o'); grid on; xlabel('N actual'); ylabel('RMSE');
            title(sprintf('RMSE - %s - RBF=%s - %s', prob.title, rbf_type, opt.epsmode));

            figure; loglog(Nactual,LINF,'-o'); grid on; xlabel('N actual'); ylabel('LINF');
            title(sprintf('LINF - %s - RBF=%s - %s', prob.title, rbf_type, opt.epsmode));

            figure; loglog(Nactual,TIME,'-o'); grid on; xlabel('N actual'); ylabel('Time');
            title(sprintf('Time - %s - RBF=%s - %s', prob.title, rbf_type, opt.epsmode));

            if strcmpi(opt.domain,'cassini')

                % ===== Affichage type surface + frontiere + points interieurs =====
                % Solution numerique :
                Uplot = Ulast;
                fig_title = sprintf('%s - %s - RBF=%s - %s',prob.title,opt.domain,rbf_type,opt.epsmode);

                % Pour afficher la solution exacte au lieu de la solution numerique,
                % remplacer les deux lignes precedentes par :
                % Uplot = prob.u(xeval(:,1),xeval(:,2));
                % fig_title = 'Profile of exact solution';

                tri = delaunay(xeval(:,1),xeval(:,2));

                % Supprimer les triangles en dehors du domaine Cassini
                % Le centre de chaque triangle doit etre calcule ligne par ligne.
                ctri_x = (xeval(tri(:,1),1) + xeval(tri(:,2),1) + xeval(tri(:,3),1))/3;
                ctri_y = (xeval(tri(:,1),2) + xeval(tri(:,2),2) + xeval(tri(:,3),2))/3;

                in_tri = inpolygon(ctri_x,ctri_y,bdpt(:,1),bdpt(:,2));
                tri = tri(in_tri(:),:);

                zmin = min(Uplot);
                zmax = max(Uplot);
                z0 = zmin - 0.25*(zmax-zmin);

                bdclose = [bdpt; bdpt(1,:)];

                figure; hold on;

                trisurf(tri,xeval(:,1),xeval(:,2),Uplot, ...
                    'EdgeColor','none','FaceColor','interp');
                shading interp;
                colormap parula;

                plot3(bdclose(:,1),bdclose(:,2),z0*ones(size(bdclose,1),1), ...
                    'b-','LineWidth',2);

                plot3(data(indI,1),data(indI,2),z0*ones(numel(indI),1), ...
                    'r.','MarkerSize',9);

                xlabel('x');
                ylabel('y');
                zlabel('v(x,y)');
                colorbar;
                grid on;
                box on;
                axis tight;
                view(135,25);
                title(fig_title);

            else
                x = linspace(0,1,opt.s_d);
                [X,Y] = meshgrid(x,x);
                figure; surf(X,Y,reshape(Ulast,opt.s_d,opt.s_d),'EdgeColor','none');
                view(3); grid on; axis square;
                title(sprintf('%s - %s - RBF=%s - %s',prob.title,opt.domain,rbf_type,opt.epsmode));
            end
        end
    end
end
end % imode
end % irep

S = groupsummary(AllResults, ...
    {'RBF','Domain','Mode','lambda','Ninput'}, ...
    {'mean','std','max'}, ...
    {'N','Eps_mean','RMSE','LINF','Time_s', ...
     'Cond_mean','Cond_max','Data_amplification'});

T = S(:,{'RBF','Domain','Mode','lambda','Ninput', ...
    'mean_N','mean_Eps_mean','mean_RMSE','std_RMSE', ...
    'mean_LINF','mean_Time_s','mean_Cond_mean', ...
    'max_Cond_max','mean_Data_amplification'});

T.Properties.VariableNames = {'RBF','Domain','Mode','lambda', ...
    'Ninput','N','Eps_mean','RMSE','RMSE_std', ...
    'LINF','Time_s','Cond_mean','Cond_max', ...
    'Data_amplification'};

disp(T);
writetable(T,'EDP_Results.csv');
% =========================================================
% Probleme
% =========================================================
function prob = get_problem(name,lambda)
switch lower(name)
    case 'poisson'
        prob.u = @(x,y) sin(pi*x).*sin(pi*y);
        prob.f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
        prob.g = @(x,y) sin(pi*x).*sin(pi*y);
        prob.title = 'Poisson';

    case 'helmholtz'
        prob.u = @(x,y) sin(pi*x).*sin(pi*y);
        prob.f = @(x,y) (2*pi^2+lambda)*sin(pi*x).*sin(pi*y);
        prob.g = @(x,y) sin(pi*x).*sin(pi*y);
        prob.title = sprintf('Helmholtz (lambda=%.2f)',lambda);

    case 'burgers_nl'
        prob.u  = @(x,y) sin(pi*x).*sin(pi*y);
        prob.ux = @(x,y) pi*cos(pi*x).*sin(pi*y);
        prob.uy = @(x,y) pi*sin(pi*x).*cos(pi*y);
        prob.f  = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y) ...
            + pi*sin(pi*x).*cos(pi*x).*sin(pi*y).^2 ...
            + pi*sin(pi*x).^2.*sin(pi*y).*cos(pi*y);
        prob.g = @(x,y) sin(pi*x).*sin(pi*y);
        prob.title = 'Nonlinear elliptic PDE';
     case 'singular_linear'
        prob.u = @(x,y) x.^2 .* cos(pi*y);
        prob.f = @(x,y) (4 - pi^2) .* cos(pi*y);
        prob.g = @(x,y) x.^2 .* cos(pi*y);
        prob.title = 'Singular linear problem';

    otherwise
        error('Probleme inconnu.');
end
end

% =========================================================
% Noeuds et grille
% =========================================================
function [data,indI,indB] = generate_nodes(N,irep)
    nb = max(20,round(sqrt(N)/2));
    t = linspace(0,1,nb)';
    B = unique([t,0*t; 1+0*t,t; flipud(t),1+0*t; 0*t,flipud(t)],'rows','stable');
    Ni = N-size(B,1);
    if Ni <= 0
        error('N doit etre superieur au nombre de points de frontiere (%d).',size(B,1));
    end
    hs = haltonset(2);
    if irep>1
        rng(1000+irep,'twister');
        hs = scramble(hs,'RR2');
    end
    X = net(hs,2*Ni);
    X = X(X(:,1)>0 & X(:,1)<1 & X(:,2)>0 & X(:,2)<1,:);
    if size(X,1) < Ni
        error('Pas assez de points interieurs.');
    end
    Xi = X(1:Ni,:);
    data = [Xi; B];
    indI = (1:size(Xi,1))';
    indB = (size(Xi,1)+1:size(data,1))';
end

function X = make_grid(s)
    [x,y] = meshgrid(linspace(0,1,s));
    X = [x(:),y(:)];
end

function [data,indI,indB,bdpt] = generate_nodes_cassini(m,nn,xshift,xzero_tol,irep)
    
    theta = 0:2*pi/nn:2*pi-2*pi/nn;
    
    r = (cos(3*theta) + sqrt(3 - sin(3*theta).^2)).^(1/3);
    
    bdpt(:,1) = r(:).*cos(theta(:)) + xshift;
    bdpt(:,2) = r(:).*sin(theta(:));
    
    a  = min(bdpt(:,1));
    b  = max(bdpt(:,1));
    c1 = min(bdpt(:,2));
    d  = max(bdpt(:,2));
    
    hx = (b-a)/m;
    hy = (d-c1)/m;
    
    if irep==1
        [M,Y] = meshgrid(a:hx:b, c1:hy:d); % grille originale
    else
        rng(1000+irep,'twister');
        dx=(rand-0.5)*0.8*hx; dy=(rand-0.5)*0.8*hy;
        [M,Y] = meshgrid((a+dx):hx:b,(c1+dy):hy:d); % decalage aleatoire reproductible
    end
    
    in = inpolygon(M(:),Y(:),bdpt(:,1),bdpt(:,2));
    intnode = [M(in), Y(in)];

    % Pour le probleme singulier, eviter les points avec x=0.
    % Si xzero_tol=0, aucun filtrage n'est applique.
    if nargin < 4
        xzero_tol = 0;
    end

    if xzero_tol > 0
        intnode = intnode(abs(intnode(:,1)) > xzero_tol,:);
        %bdpt    = bdpt(abs(bdpt(:,1)) > xzero_tol,:);
    end
    
    data = [intnode; bdpt];
    
    indI = (1:size(intnode,1))';
    indB = (size(intnode,1)+1:size(data,1))';

end

function X = make_grid_cassini(s,bdpt)

    a  = min(bdpt(:,1));
    b  = max(bdpt(:,1));
    c1 = min(bdpt(:,2));
    d  = max(bdpt(:,2));
    
    [x,y] = meshgrid(linspace(a,b,s),linspace(c1,d,s));
    
    in = inpolygon(x(:),y(:),bdpt(:,1),bdpt(:,2));
    
    X = [x(in), y(in)];

end

% =========================================================
% PU setup
% =========================================================
function PU = pu_setup(data,m_d)
    xmin = min(data(:,1)); xmax = max(data(:,1));
    ymin = min(data(:,2)); ymax = max(data(:,2));
    
    xc = linspace(xmin,xmax,m_d+1);
    yc = linspace(ymin,ymax,m_d+1);
    
    Hx = min(diff(xc));
    Hy = min(diff(yc));
    
    cx = xc(1:end-1) + Hx/2;
    cy = yc(1:end-1) + Hy/2;
    
    [cx,cy] = meshgrid(cx,cy);
    center = [cx(:),cy(:)];
    
    R = 1.5*hypot(Hx/2,Hy/2);
    
    M = size(center,1);
    patches = arrayfun(@(j) find(hypot(data(:,1)-center(j,1),data(:,2)-center(j,2))<R),1:M,'uni',0)';
    
    % Pour un domaine non rectangulaire, certains patchs de la boite peuvent etre vides.
    keep = ~cellfun(@isempty,patches);
    center = center(keep,:);
    patches = patches(keep);
    M = numel(patches);
    
    PU = struct('center',center,'R',R,'M',M,'patches',{patches});
    end
    
function PU = prepare_local_pu(data,PU)
    PU.Xloc = cell(PU.M,1);
    for j = 1:PU.M
        PU.Xloc{j} = data(PU.patches{j},:);
    end
end

function PU = attach_local_mats(PU,eps_patch,rbf_type)
    PU.Aloc = cell(PU.M,1);
    for j = 1:PU.M
        Xj = PU.Xloc{j};
        DM = distmat(Xj,Xj);
        PU.Aloc{j} = rbf_kernel_matrix(DM,eps_patch(j),rbf_type) + 1e-12*eye(size(Xj,1));
        %PU.Aloc{j} = rbf_kernel_matrix(DM,eps_patch(j),rbf_type);

    end
end

% =========================================================
% Assemblage
% =========================================================
function [Axx,Ayy] = build_pu_laplacian(data,PU,eps_patch,rbf_type)
    
    N = size(data,1);
    Axx = sparse(N,N);
    Ayy = sparse(N,N);
    
    for i = 1:N
        x = data(i,:);
        act = active_patches(x,PU);
        [w,wx,wy,wxx,wyy] = weight_wendland_c2(x,PU.center(act,:),PU.R);
    
        rx = sparse(1,N);
        ry = sparse(1,N);
    
        for k = 1:numel(act)
            j = act(k);
            idx = PU.patches{j};
            [psi,psix,psiy,psixx,psiyy] = rbf_ops_2d(eps_patch(j),PU.Xloc{j},x,PU.Aloc{j},rbf_type);
    
            rx(idx) = rx(idx) + wxx(k)*psi + 2*wx(k)*psix + w(k)*psixx;
            ry(idx) = ry(idx) + wyy(k)*psi + 2*wy(k)*psiy + w(k)*psiyy;
        end
    
        Axx(i,:) = rx;
        Ayy(i,:) = ry;
    end
    end
    
function [Ax,Ay,Axx,Ayy] = build_pu_operators(data,PU,eps_patch,rbf_type)
    N = size(data,1);
    Ax  = sparse(N,N);
    Ay  = sparse(N,N);
    Axx = sparse(N,N);
    Ayy = sparse(N,N);
    
    for i = 1:N
        x = data(i,:);
        act = active_patches(x,PU);
        [w,wx,wy,wxx,wyy] = weight_wendland_c2(x,PU.center(act,:),PU.R);
    
        rx  = sparse(1,N);
        ry  = sparse(1,N);
        rxx = sparse(1,N);
        ryy = sparse(1,N);
    
        for k = 1:numel(act)
            j = act(k);
            idx = PU.patches{j};
            [psi,psix,psiy,psixx,psiyy] = rbf_ops_2d(eps_patch(j),PU.Xloc{j},x,PU.Aloc{j},rbf_type);
    
            rx(idx)  = rx(idx)  + wx(k)*psi  + w(k)*psix;
            ry(idx)  = ry(idx)  + wy(k)*psi  + w(k)*psiy;
            rxx(idx) = rxx(idx) + wxx(k)*psi + 2*wx(k)*psix + w(k)*psixx;
            ryy(idx) = ryy(idx) + wyy(k)*psi + 2*wy(k)*psiy + w(k)*psiyy;
        end
    
        Ax(i,:)  = rx;
        Ay(i,:)  = ry;
        Axx(i,:) = rxx;
        Ayy(i,:) = ryy;
    end
end

function U = eval_pu_solution(xeval,data,u,PU,eps_patch,rbf_type)
    U = zeros(size(xeval,1),1);
    for i = 1:size(xeval,1)
        x = xeval(i,:);
        act = active_patches(x,PU);
        [w,~,~,~,~] = weight_wendland_c2(x,PU.center(act,:),PU.R);
    
        for k = 1:numel(act)
            j = act(k);
            idx = PU.patches{j};
            psi = rbf_ops_2d(eps_patch(j),PU.Xloc{j},x,PU.Aloc{j},rbf_type);
            U(i) = U(i) + w(k)*(psi*u(idx));
        end
    end
end

function act = active_patches(x,PU)
    act = find(hypot(x(1)-PU.center(:,1),x(2)-PU.center(:,2)) < PU.R);
    
    % Petite securite numerique si un point tombe exactement hors recouvrement.
    if isempty(act)
        [~,act] = min(hypot(x(1)-PU.center(:,1),x(2)-PU.center(:,2)));
    end
end

% =========================================================
% LOOCV EDP de type c_opt, patch par patch
% =========================================================
function copt = optimize_c(obj,base_obj,bounds,tolx,mode)
    options=optimset('TolX',tolx);
    if strcmpi(mode,'loocv_ga')
        if exist('ga','file')~=2
            error('Mode loocv_ga : Global Optimization Toolbox requis.');
        end
        go=optimoptions('ga','PopulationSize',16,'MaxGenerations',12,'Display','off');
        copt=ga(@(z)obj(z(1)),1,[],[],[],[],bounds(1),bounds(2),[],go);
        return
    end
    if ismember(lower(mode),{'loocv','econm_c'})
        copt=fminbnd(obj,bounds(1),bounds(2),options);
        return
    end
    % Comparateur: recherche deterministe sur grille+raffinement local.
    v=linspace(bounds(1),bounds(2),17);
    cost=arrayfun(obj,v);
    cost(~isfinite(cost) | cost>realmax/10)=Inf;
    if all(isinf(cost))
        % Aucun candidat admissible: conserver le meilleur LOOCV NON contraint.
        % Les conditionnements non admissibles sont explicites dans les tableaux.
        cost=arrayfun(base_obj,v);
        cost(~isfinite(cost) | cost>realmax/10)=Inf;
        obj=base_obj;
    end
    [best,k]=min(cost);
    if ~isfinite(best), copt=mean(bounds); return; end
    a=v(max(1,k-1)); b=v(min(numel(v),k+1));
    if a<b
        [t,ft]=fminbnd(obj,a,b,options);
        if isfinite(ft) && ft<best, copt=t; else, copt=v(k); end
    else
        copt=v(k);
    end
end

function ep = choose_eps_pde_patch(X,isB,prob,problem,lambda,c_bounds,tolx,rbf_type,mode,cond_limit)

    if size(X,1) < 5
        ep = parameter_to_epsilon(mean(c_bounds),rbf_type);
        return;
    end

    obj = @(c) loocv_cost_pde_patch(c,X,isB,prob,problem,lambda,rbf_type,mode,cond_limit);
    base_obj = @(c) loocv_cost_pde_patch(c,X,isB,prob,problem,lambda,rbf_type,'loocv',cond_limit);
    copt = optimize_c(obj,base_obj,c_bounds,tolx,mode);
    ep = parameter_to_epsilon(copt,rbf_type);
end

function s = loocv_cost_pde_patch(c,X,isB,prob,problem,lambda,rbf_type,mode,cond_limit)
    ep = parameter_to_epsilon(c,rbf_type);
    n = size(X,1);
    C = zeros(n,n);
    rhs = zeros(n,1);
    I = eye(n);

    for q = 1:n
        xq = X(q,:);

        if isB(q)
            % Condition de Dirichlet : u=g sur la frontiere.
            psi = rbf_ops_2d(ep,X,xq,I,rbf_type);
            C(q,:) = psi;
            rhs(q) = prob.g(xq(1),xq(2));
        else
            [psi,psix,~,psixx,psiyy] = rbf_ops_2d(ep,X,xq,I,rbf_type);

            switch lower(problem)
                case 'poisson'
                    % -Delta u = f
                    C(q,:) = -(psixx + psiyy);

                case 'helmholtz'
                    % -Delta u + lambda*u = f
                    C(q,:) = -(psixx + psiyy) + lambda*psi;

                case 'singular_linear'
                    % u_xx + (1/x^2)u_yy + (1/x)u_x = f
                    if abs(xq(1)) < 1e-12
                        s = realmax;
                        return;
                    end
                    C(q,:) = psixx + (1/xq(1)^2)*psiyy + (1/xq(1))*psix;

                otherwise
                    % Le code c_opt fourni est base sur un operateur lineaire.
                    s = realmax;
                    return;
            end

            rhs(q) = prob.f(xq(1),xq(2));
        end
    end

    if strcmpi(mode,'econm_c')
        s = econm_cost_patch(C,X,ep,rbf_type);
        return;
    end
    if strcmpi(mode,'loocv_cond')
        Kr=rbf_kernel_matrix(distmat(X,X),ep,rbf_type);
        if rcond(Kr)<1/cond_limit || rcond(C)<1/cond_limit
            s=Inf; return
        end
    end

    % Meme petite regularisation que dans les matrices locales du main.
    C = C + 1e-12*I;

    % Critere LOOCV de meme forme que dans c_opt :
    % errorvector = coefficients ./ diag(inv(C)), sans former inv(C).
    coeff = C\rhs;
    dinv = diag(C\I);

    if any(~isfinite(coeff)) || any(~isfinite(dinv)) || any(abs(dinv) < 1e-14)
        s = realmax;
        return;
    end

    e = coeff./dinv;
    s = norm(e,2);

    if ~isfinite(s)
        s = realmax;
    end
end

% =========================================================
% LOOCV Burgers : meme construction RBF brute que les EDP lineaires
% =========================================================
function ep = choose_eps_burgers_patch(X,isB,prob,c_bounds,tolx,rbf_type,mode,cond_limit,varargin)
    % Premier appel : varargin vide -> ubar=uxbar=uybar=0 -> -Delta.
    % Deuxieme appel : varargin = {ubar,uxbar,uybar} -> operateur linearise.
    if size(X,1) < 5
        ep = parameter_to_epsilon(mean(c_bounds),rbf_type);
        return;
    end

    ubar = [];
    uxbar = [];
    uybar = [];
    if numel(varargin) >= 3
        ubar = varargin{1};
        uxbar = varargin{2};
        uybar = varargin{3};
    end

    obj = @(c) loocv_cost_burgers_patch( ...
        c,X,isB,prob,rbf_type,ubar,uxbar,uybar,mode,cond_limit);
    base_obj = @(c) loocv_cost_burgers_patch( ...
        c,X,isB,prob,rbf_type,ubar,uxbar,uybar,'loocv',cond_limit);
    copt = optimize_c(obj,base_obj,c_bounds,tolx,mode);
    ep = parameter_to_epsilon(copt,rbf_type);
end

function s = loocv_cost_burgers_patch(c,X,isB,prob,rbf_type,ubar,uxbar,uybar,mode,cond_limit)
    ep = parameter_to_epsilon(c,rbf_type);
    n = size(X,1);
    C = zeros(n,n);
    rhs = zeros(n,1);
    I = eye(n);

    for q = 1:n
        xq = X(q,:);

        if isB(q)
            % Dirichlet : u=g. Utiliser la RBF brute, comme dans la
            % premiere LOOCV lineaire validee (A=I dans rbf_ops_2d).
            psi = rbf_ops_2d(ep,X,xq,I,rbf_type);
            C(q,:) = psi;
            rhs(q) = prob.g(xq(1),xq(2));
        else
            [psi,psix,psiy,psixx,psiyy] = ...
                rbf_ops_2d(ep,X,xq,I,rbf_type);

            if isempty(ubar)
                % Premiere selection : ubar=0 -> J_0(v)=-Delta v.
                C(q,:) = -(psixx + psiyy);
                rhs(q) = prob.f(xq(1),xq(2));
            else
                % Linearisation de Newton autour de la solution numerique ubar :
                % J_ubar(v) = -Delta v +(uxbar+uybar)v + ubar(v_x+v_y).
                C(q,:) = -(psixx + psiyy) ...
                    + (uxbar(q)+uybar(q))*psi ...
                    + ubar(q)*(psix+psiy);

                % J_ubar(u_new) = f + ubar(uxbar+uybar).
                rhs(q) = prob.f(xq(1),xq(2)) ...
                    + ubar(q)*(uxbar(q)+uybar(q));
            end
        end
    end

    if strcmpi(mode,'econm_c')
        s = econm_cost_patch(C,X,ep,rbf_type);
        return;
    end
    if strcmpi(mode,'loocv_cond')
        Kr=rbf_kernel_matrix(distmat(X,X),ep,rbf_type);
        if rcond(Kr)<1/cond_limit || rcond(C)<1/cond_limit
            s=Inf; return
        end
    end

    C = C + 1e-12*I;
    coeff = C\rhs;
    dinv = diag(C\I);

    if any(~isfinite(coeff)) || any(~isfinite(dinv)) || any(abs(dinv) < 1e-14)
        s = realmax;
        return;
    end

    e = coeff./dinv;
    s = norm(e,2);
    if ~isfinite(s)
        s = realmax;
    end
end

% =========================================================
% Critere econm_c 
function score = econm_cost_patch(A,X,ep,rbf_type)
    [~,kk] = min(sum((X-mean(X,1)).^2,2));
    D = sqrt(sum((X-X(kk,:)).^2,2));
    rhs = rbf_kernel_matrix(D,ep,rbf_type);
    if any(~isfinite(A(:))) || any(~isfinite(rhs(:))) || norm(rhs)==0
        score = realmax;
        return;
    end
    x = lsqminnorm(A,rhs);
    S = svd(A);
    S_min = S(end);
    score = S_min*norm(x)/norm(rhs); % Equivalent aux deux lignes de CostR
    if ~isfinite(score)
        score = realmax;
    end
end

% =========================================================
function ep = parameter_to_epsilon(c,rbf_type)
    ep = c;
end

% =========================================================
% RBF kernels et derivees
% =========================================================
function check_rbf_type(rbf_type)
    if ~ismember(lower(rbf_type),{'sinc','mq','gaussian'})
        error('Noyau RBF inconnu. Utiliser ''sinc'', ''mq'' ou ''gaussian''.');
    end
end

function A = rbf_kernel_matrix(DM,ep,rbf_type)
    switch lower(rbf_type)
        case 'sinc'
            % MATLAB sinc(z)=sin(pi*z)/(pi*z)
            A = sinc(ep*DM);
    
        case 'mq'
            % Multiquadric : phi(r)=sqrt(1+(epsilon*r)^2)
            A = sqrt(1 + (ep*DM).^2);

        case 'gaussian'
            A = exp(-(ep*DM).^2);

        otherwise
            error('Noyau RBF inconnu.');
    end
end

function [psi,psix,psiy,psixx,psiyy] = rbf_ops_2d(ep,xc,xe,A,rbf_type)
% Calcule les cardinales locales et leurs derivees pour :
% - sinc : phi(r)=sin(pi*epsilon*r)/(pi*epsilon*r)
% - mq   : phi(r)=sqrt(1+(epsilon*r)^2)

    n = size(xc,1);
    dx = xe(1)-xc(:,1);
    dy = xe(2)-xc(:,2);
    r = sqrt(dx.^2+dy.^2);
    nz = r > 1e-14;
    
    phi   = zeros(n,1);
    phir  = zeros(n,1);
    phirr = zeros(n,1);
    
    switch lower(rbf_type)
        case 'sinc'
            a = pi*ep;
            phi(:) = 1;
            phi(nz) = sin(a*r(nz))./(a*r(nz));
    
            phir(nz) = cos(a*r(nz))./r(nz) - sin(a*r(nz))./(a*r(nz).^2);
    
            phirr(nz) = -a*sin(a*r(nz))./r(nz) ...
                - 2*cos(a*r(nz))./r(nz).^2 ...
                + 2*sin(a*r(nz))./(a*r(nz).^3);
            phirr(~nz) = -a^2/3;
    
        case 'mq'
            % phi(r)=sqrt(1+(ep*r)^2)
            q = 1 + (ep*r).^2;
            phi = sqrt(q);
    
            phir(nz) = ep^2*r(nz)./sqrt(q(nz));
            phirr(nz) = ep^2./sqrt(q(nz)) - ep^4*r(nz).^2./(q(nz).^(3/2));
            phirr(~nz) = ep^2;

        case 'gaussian'
            phi = exp(-(ep*r).^2);
            phir(nz) = -2*ep^2*r(nz).*phi(nz);
            phirr = (4*ep^4*r.^2 - 2*ep^2).*phi;

        otherwise
            error('Noyau RBF inconnu.');
    end
    
    phix = zeros(n,1);
    phiy = zeros(n,1);
    phixx = zeros(n,1);
    phiyy = zeros(n,1);
    
    phix(nz) = phir(nz).*dx(nz)./r(nz);
    phiy(nz) = phir(nz).*dy(nz)./r(nz);
    
    phixx(nz) = phirr(nz).*(dx(nz).^2./r(nz).^2) ...
        + phir(nz).*(1./r(nz) - dx(nz).^2./r(nz).^3);
    phiyy(nz) = phirr(nz).*(dy(nz).^2./r(nz).^2) ...
        + phir(nz).*(1./r(nz) - dy(nz).^2./r(nz).^3);
    
    % Au centre r=0 : phi_xx(0)=phi_yy(0)=phi''(0)
    phixx(~nz) = phirr(~nz);
    phiyy(~nz) = phirr(~nz);
    
    psi = phi'/A;
    
    if nargout > 1
        psix  = phix'/A;
        psiy  = phiy'/A;
        psixx = phixx'/A;
        psiyy = phiyy'/A;
    end
end

% =========================================================
% Outils
% =========================================================
function DM = distmat(X,Y)
    DM = sqrt(max(sum(X.^2,2) + sum(Y.^2,2)' - 2*(X*Y'),0));
end

function [w,wx,wy,wxx,wyy] = weight_wendland_c2(x,cen,R)
    dx = x(1)-cen(:,1);
    dy = x(2)-cen(:,2);
    r = hypot(dx./R,dy./R);
    
    phi = (4*r+1).*(1-r).^4;
    phix = 20*dx.*(r-1).^3/R^2;
    phiy = 20*dy.*(r-1).^3/R^2;
    
    phixx = zeros(size(r));
    phiyy = phixx;
    iz = r==0;
    nz = ~iz;
    
    phixx(iz) = -20/R^2;
    phiyy(iz) = -20/R^2;
    
    phixx(nz) = 20*(r(nz)-1).^2.*(-r(nz)+4*(dx(nz).^2/R^2)+(dy(nz).^2/R^2))./(R^2*r(nz));
    phiyy(nz) = 20*(r(nz)-1).^2.*(-r(nz)+(dx(nz).^2/R^2)+4*(dy(nz).^2/R^2))./(R^2*r(nz));
    
    s = sum(phi);
    sx = sum(phix);
    sy = sum(phiy);
    sxx = sum(phixx);
    syy = sum(phiyy);
    
    w = phi/s;
    wx = phix/s - phi*sx/s^2;
    wy = phiy/s - phi*sy/s^2;
    wxx = phixx/s - 2*phix*sx/s^2 + phi*(2*sx^2/s^3 - sxx/s^2);
    wyy = phiyy/s - 2*phiy*sy/s^2 + phi*(2*sy^2/s^3 - syy/s^2);
end
