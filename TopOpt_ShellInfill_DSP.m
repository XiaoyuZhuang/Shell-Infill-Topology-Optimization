%% TopOpt_ShellInfill_DSP
% -------------------------------------------------------------------------
% Description:
%   Topology Optimization for Shell-Infill Structures using the DSP 
%   (Design Space Projection) and Erosion-based method.
%
% Key Features:
%   1. Double-Filter approach: Coarse structure (DSP) + Erosion (Shell).
%   2. Gradient Descent with Adaptive Step (GDA) optimizer.
%   3. Special Boundary Condition: Open/Symmetric boundary on the left edge.
% -------------------------------------------------------------------------

function TopOpt_ShellInfill_DSP
    clc; clear; close all;
    
    %% === 1. CONFIGURATION & PARAMETERS ===
    % Domain and Mesh
    nelx = 120;             % Number of elements in X
    nely = 60;              % Number of elements in Y
    volfrac = 0.5;          % Target volume fraction
    
    % Optimization Settings
    penal = 3;              % Penalization factor (SIMP)
    step_size = 0.1;        % Gradient descent step size
    move_limit = 0.05;      % Move limit per iteration
    max_loop = 3000;         % Maximum iterations (safety break)
    tol_ch = 0.005;         % Convergence tolerance
    
    % Projection & Filter Parameters
    R_dsp = 6.0;            % Filter radius for DSP (Base structure size)
    R_ero = 3.0;            % Erosion radius (Shell thickness control)
    eta_dsp = 0.50;         % Projection threshold for base structure
    eta_ero = 0.80;         % Projection threshold for erosion
    beta = 1;               % Initial projection steepness
    beta_max = 64;          % Maximum projection steepness
    beta_inc_interval = 100;% Interval to double beta
    
    % Material Properties (Normalized)
    E_shell = 1.0;          % Young's modulus: Shell
    E_infill = 0.2;         % Young's modulus: Infill
    E_min = 1e-9;           % Young's modulus: Void
    nu = 0.3;               % Poisson's ratio
    
    %% === 2. FINITE ELEMENT PRE-CALCULATION ===
    % Compute element stiffness matrix (KE)
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
    
    % Index maps for assembly
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
    iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
    
    % Boundary Conditions and Loads (Cantilever with symmetry on Left)
    % Fixed nodes: Left edge restricted, logic defined by user
    fixeddofs = union([1:2:2*(nely+1)], [2*(nelx+1)*(nely+1)]);
    alldofs = 1:2*(nely+1)*(nelx+1);
    freedofs = setdiff(alldofs, fixeddofs);
    
    % Load: Point load at bottom-right corner
    F = sparse(2, 1, -1, 2*(nely+1)*(nelx+1), 1);
    U = zeros(2*(nely+1)*(nelx+1), 1);
    
    %% === 3. FILTER INITIALIZATION ===
    % H1: DSP Filter (Base Structure)
    [dy, dx] = meshgrid(-ceil(R_dsp)+1 : ceil(R_dsp)-1, -ceil(R_dsp)+1 : ceil(R_dsp)-1);
    h1 = max(0, R_dsp - sqrt(dx.^2 + dy.^2));
    Hs_dsp = conv2(ones(nely, nelx), h1, 'same');
    
    % H2: Erosion Filter (Constant kernel for geometric erosion)
    pad_size = ceil(R_ero); 
    [dy, dx] = meshgrid(-pad_size:pad_size, -pad_size:pad_size);
    h2 = max(0, R_ero - sqrt(dx.^2 + dy.^2));
    Hs_ero_scalar = sum(h2(:)); % Scalar normalization for uniform erosion
    
    %% === 4. OPTIMIZATION LOOP ===
    x = repmat(volfrac, nely, nelx); 
    loop = 0; 
    change = 1;
    
    % Initialize Figures
    figure(1); set(gcf, 'Name', 'Convergence History', 'Position', [100, 100, 600, 400]);
    figure(2); set(gcf, 'Name', 'Process Visualization', 'Position', [750, 100, 1000, 600]);
    
    while (change > tol_ch || beta < beta_max) && loop < max_loop
        loop = loop + 1;
        
        %% --- A. Forward Projection (Physical Fields) ---
        % 1. DSP Projection: Design Variable (x) -> Base Topology (mu)
        xTilde = conv2(x, h1, 'same') ./ Hs_dsp;
        [xBar, dxBar_dxT] = heaviside_proj(xTilde, beta, eta_dsp);
        xHat = conv2(xBar, h1, 'same') ./ Hs_dsp;
        [mu, dmu_dxH] = heaviside_proj(xHat, beta, eta_dsp);
        
        % 2. Erosion Projection: Base Topology (mu) -> Infill Core (phi)
        % Padding: Handle left boundary (Symmetry/No-shell condition)
        mu_aug = [repmat(mu(:,1), 1, pad_size), mu]; 
        
        % Convolution with H2
        muTilde_aug = conv2(mu_aug, h2, 'same') ./ Hs_ero_scalar;
        muTilde = muTilde_aug(:, pad_size+1 : end); 
        
        % Projection to get Core
        [phi, dphi_dmuT] = heaviside_proj(muTilde, beta, eta_ero);
        
        % 3. Material Distribution
        rho_shell = max(0, mu - phi); % Shell region
        rho_infill = phi;             % Infill region
        
        % SIMP Interpolation
        E_elem = E_min + rho_shell.^penal * E_shell + rho_infill.^penal * E_infill;
        
        %% --- B. Visualization ---
        if mod(loop, 1) == 0
            plot_process(loop, beta, mu, muTilde, phi, rho_shell, rho_infill, x, eta_ero);
        end
        
        %% --- C. Finite Element Analysis ---
        sK = reshape(KE(:)*E_elem(:)', 64*nelx*nely, 1);
        K = sparse(iK, jK, sK); K = (K+K')/2;
        U(freedofs) = K(freedofs,freedofs) \ F(freedofs);
        
        % Objective Function (Compliance)
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2), nely, nelx);
        c_obj = sum(sum(E_elem.*ce));
        
        %% --- D. Sensitivity Analysis ---
        dc_dE = -ce;
        
        % Derivatives w.r.t densities
        dE_drho_shell = penal * rho_shell.^(penal-1) * E_shell;
        dE_drho_infill = penal * rho_infill.^(penal-1) * E_infill;
        
        % Derivatives w.r.t physical fields (mu, phi)
        dE_dmu_part = dE_drho_shell; 
        dE_dphi = -dE_drho_shell + dE_drho_infill;
        
        % Chain rule for Erosion path (phi -> muTilde -> mu)
        sens_phi = dc_dE .* dE_dphi;
        sens_muTilde = sens_phi .* dphi_dmuT;
        
        % Handle padding sensitivity (Adjoint operation of padding)
        sens_muTilde_aug = zeros(nely, nelx + pad_size);
        sens_muTilde_aug(:, pad_size+1 : end) = sens_muTilde;
        sens_mu_aug = conv2(sens_muTilde_aug ./ Hs_ero_scalar, h2, 'same');
        
        sens_mu_erosion = sens_mu_aug(:, pad_size+1 : end); 
        sens_mu_erosion_left = sum(sens_mu_aug(:, 1:pad_size), 2);
        sens_mu_erosion(:, 1) = sens_mu_erosion(:, 1) + sens_mu_erosion_left;
        
        % Combine sensitivities on mu
        sens_mu_direct = dc_dE .* dE_dmu_part;
        dc_dmu = sens_mu_direct + sens_mu_erosion;
        
        % Chain rule for DSP path (mu -> xHat -> xBar -> xTilde -> x)
        sens_xHat = dc_dmu .* dmu_dxH;        
        sens_xBar = conv2(sens_xHat ./ Hs_dsp, h1, 'same'); 
        sens_xTilde = sens_xBar .* dxBar_dxT; 
        dc_dx = conv2(sens_xTilde ./ Hs_dsp, h1, 'same');   
        
        % Volume Sensitivity
        dv_dmu = ones(nely,nelx) / (nelx*nely);
        dv_dxHat = dv_dmu .* dmu_dxH;
        dv_dxBar = conv2(dv_dxHat ./ Hs_dsp, h1, 'same');
        dv_dxTilde = dv_dxBar .* dxBar_dxT;
        dv_dx = conv2(dv_dxTilde ./ Hs_dsp, h1, 'same');
        
        current_vol = mean(mu(:));
        
        %% --- E. Design Update (GDA with Volume Constraint) ---
        % Normalize gradient
        scaling = max(abs(dc_dx(:)));
        if scaling > 0, dc_dx_norm = dc_dx / scaling; else, dc_dx_norm = dc_dx; end
        
        % Bisection method to find Lagrange multiplier (lambda)
        l1 = 0; l2 = 1e8; 
        x_new = x;
        
        while (l2-l1)/(l1+l2) > 1e-4
            lmid = 0.5*(l2+l1);
            
            % Update direction
            change_step = -step_size * (dc_dx + lmid * dv_dx) ./ max(1e-10, max(abs(dc_dx(:) + lmid * dv_dx(:))));
            x_trial = x + change_step;
            
            % Apply Move Limit and Bounds
            x_trial = max(0, min(1, x_trial));
            x_trial = max(x - move_limit, min(x + move_limit, x_trial));
            
            % Predict Volume
            xT_new = conv2(x_trial, h1, 'same') ./ Hs_dsp;
            xB_new = heaviside_proj(xT_new, beta, eta_dsp);
            xH_new = conv2(xB_new, h1, 'same') ./ Hs_dsp;
            mu_new = heaviside_proj(xH_new, beta, eta_dsp);
            
            if mean(mu_new(:)) > volfrac
                l1 = lmid; 
            else
                l2 = lmid; 
            end
        end
        x_new = x_trial;
        change = max(abs(x_new(:)-x(:)));
        x = x_new;
        
        %% --- F. Continuation Scheme ---
        if beta < beta_max && (mod(loop, beta_inc_interval) == 0 || change < tol_ch)
            beta = min(beta_max, beta * 2);
            fprintf('>>> Beta updated to %d\n', beta);
        end
        
        fprintf(' It.:%4i Obj.:%10.4f Vol.:%6.3f ch.:%6.3f Lam:%6.2e\n', ...
            loop, c_obj, current_vol, change, lmid);
        
        % Update Status Figure
        set(0, 'CurrentFigure', 1);
        viz = rho_shell * 1.0 + rho_infill * 0.5;
        imagesc(-viz); axis equal tight off; colormap(gray);
        title(['Iteration: ' num2str(loop) ', Beta: ' num2str(beta)]); drawnow;
    end
end

%% Helper Functions
function [rho, drho_dt] = heaviside_proj(t, beta, eta)
    num = tanh(beta * eta) + tanh(beta * (t - eta));
    den = tanh(beta * eta) + tanh(beta * (1 - eta));
    rho = num ./ den;
    dnum = beta * (1 - tanh(beta * (t - eta)).^2);
    drho_dt = dnum ./ den;
end

function plot_process(loop, beta, mu, muTilde, phi, rho_shell, rho_infill, x, eta_ero)
    set(0, 'CurrentFigure', 2);
    
    subplot(2,3,1); 
    imagesc(-mu); axis equal tight off; colormap(gray); 
    title('1. Base Structure (\mu)');
    
    subplot(2,3,2); 
    imagesc(-muTilde); axis equal tight off; 
    title('2. Eroded Gradient (\mu_{tilde})');
    hold on; contour(muTilde, [eta_ero eta_ero], 'r', 'LineWidth', 1); hold off;
    
    subplot(2,3,3); 
    imagesc(-phi); axis equal tight off; 
    title('3. Core Phase (\phi)');
    
    subplot(2,3,4); 
    imagesc(-rho_shell); axis equal tight off; 
    title('4. Shell Phase');
    
    subplot(2,3,5); 
    imagesc(-rho_infill); axis equal tight off; 
    title('5. Infill Phase');
    
    subplot(2,3,6); 
    imagesc(-x); axis equal tight off; 
    title('6. Design Variable (x)');
    
    sgtitle(['Iteration: ' num2str(loop) ' | Beta: ' num2str(beta)]);
    drawnow;
end
