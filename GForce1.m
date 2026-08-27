%% GForce1
function [kt, kn, kr, s_ang_final, tau1_final, Vchip_final, mu_final, f_ang_final, p0, px, h_cr, signz, tau_f, l1, l2, l3, time] = GForce1(gm, gs, chip, inc, vc, Tr, Tins, material)
    % This function calculates the forces involved in the operation in
    % microscale

    t = tic;
    
    %% Stage 1, cutting force coefficients
    %syms btn
    diskrete = 50;
    ktc = cell(size(gm, 1), size(gm, 2));
    kt = zeros(size(gm, 1), size(gm, 2));
    kfc = cell(size(gm, 1), size(gm, 2));
    kn = zeros(size(gm, 1), size(gm, 2));
    krc = cell(size(gm, 1), size(gm, 2));
    kr = zeros(size(gm, 1), size(gm, 2));
    gamma1 = zeros(size(gm, 1), size(gm, 2), 20);
    % signGamma1 = zeros(size(gm, 1), size(gm, 2), 20);
    gamma_dot = zeros(size(gm, 1), size(gm, 2), 20);
    tau = zeros(size(gm, 1), size(gm, 2), 20);
    eta_s = zeros(size(gm, 1), size(gm, 2));
    Vchip = zeros(size(gm, 1), size(gm, 2));
    mu = zeros(size(gm, 1), size(gm, 2));
    tau1 = zeros(size(gm, 1), size(gm, 2), 20);
    beta_n = zeros(size(gm, 1), size(gm, 2), 20);
    s_ang_final = zeros(size(gm, 1), size(gm, 2));
    tau1_final = zeros(size(gm, 1), size(gm, 2));
    Vchip_final = zeros(size(gm, 1), size(gm, 2));
    mu_final = zeros(size(gm, 1), size(gm, 2));
    f_ang_final = zeros(size(gm, 1), size(gm, 2));
    Lcf = zeros(size(gm, 1), size(gm, 2));
    % tau1Overp0 = zeros(1, diskrete);
    % signTau1Overp0 = zeros(1, diskrete);
    tau1Overp0f = zeros(size(gm, 1), size(gm, 2));
    p0 = zeros(size(gm, 1), size(gm, 2));
    p00 = zeros(size(gm, 1), size(gm, 2));
    
    h_cr = zeros(size(gm, 1), size(gm, 2));
    l1 = zeros(size(gm, 1), size(gm, 2));
    l2 = zeros(size(gm, 1), size(gm, 2));
    l3 = zeros(size(gm, 1), size(gm, 2));
    lf = zeros(size(gm, 1), size(gm, 2));
    x_f = cell(size(gm, 1), size(gm, 2));
    a = zeros(size(gm, 1), size(gm, 2));
    b = zeros(size(gm, 1), size(gm, 2));
    c = zeros(size(gm, 1), size(gm, 2));
    px = cell(size(gm, 1), size(gm, 2));
    lst = cell(size(gm, 1), size(gm, 2));
    signz = cell(size(gm, 1), size(gm, 2));
    interindi = cell(size(gm, 1), size(gm, 2));
    tau_f = cell(size(gm, 1), size(gm, 2));
    

    chip = cellfun(@max, chip, 'UniformOutput', false);
    chip = cellfun(@(c) sum([0, c]), chip);
    % vc cutting speed m/s%
    vc = vc ./ 1e3;
    offset = 5;
    for d = 1 : 1 : size(gm, 2)
        for r = inc{1, d}
            for i = 1 : 1 : 17
                phi_n = i + offset;                
                zeta = 2;
                alpha_n = gs(r, d, 1);
                o_ang = gs(r, d, 7);
                % Stabler's rule
                eta_c = o_ang;
                
                T = Tins;

                gamma1(r, d, i) = cotd(phi_n) + tand(phi_n - alpha_n);
                signGamma1 = sign(gamma1);
                h1 = .00001; % [m]

                % Strain rate, gamma_dot, equation can lead to complexity in all other equations, if the sum of test shear angle and rake angle exceed 90 
                % abs() used for cosd to prevent negative number for strain rate and complexity in tau and etc
                gamma_dot(r, d, i) = (vc(r, d) .* cosd(alpha_n)) ./ (h1 .* abs(cosd(phi_n - alpha_n)));
                gamma_dot_ref = .001;
                % Johnson-Cook
                tau(r, d, i) = 1 ./ sqrt(3) .* (material.A + material.B .* signGamma1(r, d, i) .* (abs(gamma1(r, d, i)) ./ sqrt(3)) .^ material.n) .* (1 + material.C .* log(gamma_dot(r, d, i) ./ gamma_dot_ref)) .* (1 - ((T - Tr) ./ ((material.Tm - 273.15) - Tr)) .^ material.m);
            
                eta_s(r, d) = atand((tand(o_ang) .* cosd(phi_n - alpha_n) - tand(o_ang) .* sind(phi_n)) ./ cosd(alpha_n));
                Vchip(r, d, i) = (vc(r, d) .* sind(phi_n)) ./ cosd(phi_n - alpha_n);
                mu(r, d) = .1484 + .0061 .* Vchip(r, d);
                tau1(r, d, i) = tau(r, d, i) + (1e-6 * material.rho .* ((vc(r, d) .* sind(phi_n) .* cosd(o_ang)) .^ 2) .* gamma1(r, d, i));
                
                % btn is an abbreviation for  beta_n
                btn = linspace(-10, 10, diskrete);
                Lc = (chip(r, d) .* (zeta + 2) .* sind(phi_n + btn - alpha_n)) ./ (2 .* sind(phi_n) .* cosd(btn) .* cosd(eta_c));
                tau1Overp0 = (Lc .* sind(phi_n) .* cosd(eta_c) .* cosd(phi_n + btn - alpha_n)) ./ (chip(r, d) .* (zeta + 1) .* cosd(eta_s(r, d)) .* cosd(btn));
                signTau1Overp0Overmu = sign(tau1Overp0 ./ mu(r, d));
                mu_a = tau1Overp0 .* (1 + zeta .* (1 - (signTau1Overp0Overmu .* (abs(tau1Overp0 ./ mu(r, d)) .^ (1 / zeta) ) ) ) );
                f = mu_a .* cosd(eta_c) - btn;
                efed = f;
                %btn = btn(f == 0);

                % xaxis = zeros(size(efed));
                % fplot(f, 'Color', [rand, rand, rand])
                % plot(btn, efed, ':', 'LineWidth', 2)
                % text(i, i * .1, strcat('s', int2str(i)))
                % hold on
                % grid minor
                %[beta_n(r, d, i), ~] = fzero(f, 2);
                %[IntPx, IntPz] = intersections(btn, efed, btn, xaxis);
                beta_n(r, d, i) = interp1(efed, btn, 0, 'linear');
                %beta_n(r, d, i) = IntPx;
                %f_ang(r, d, i) = tand(lambda_a(r, d, i)) * cosd(o_ang);
                f_ang = beta_n;

                denominator = sqrt(cosd(phi_n + beta_n(r, d, i) - alpha_n) .^ 2 + tand(o_ang) .^ 2 .* sind(beta_n(r, d, i)) .^ 2);
                ktc{r, d} = [ktc{r, d}, ((tau1(r, d, i) ./ sind(phi_n)) .* cosd(beta_n(r, d, i) - alpha_n) + tand(o_ang) .^ 2 .* sind(beta_n(r, d, i))) ./ denominator];
                kfc{r, d} = [kfc{r, d}, ((tau1(r, d, i) ./ sind(phi_n) .* cosd(o_ang)) .* sind(beta_n(r, d, i) - alpha_n)) ./ denominator];
                krc{r, d} = [krc{r, d}, ((tau1(r, d, i) ./ sind(phi_n)) .* (cosd(beta_n(r, d, i) - alpha_n) .* tand(o_ang) - tand(o_ang) .* sind(beta_n(r, d, i)))) ./ denominator];

            end

            ktc{r, d} = ktc{r, d};
            kfc{r, d} = kfc{r, d};
            krc{r, d} = krc{r, d};
            
            ktc{r, d}(ktc{r, d} < 0)  = mean(ktc{r, d}(ktc{r, d} > 5e3));
            kfc{r, d}(kfc{r, d} < 0)  = mean(kfc{r, d}(kfc{r, d} > 5e1));
            krc{r, d}(krc{r, d} < 0)  = mean(krc{r, d}(krc{r, d} > 5e3));

            [~, index] = min(ktc{r, d}, [], 'linear');
            s_ang_final(r, d) = index + offset;
            tempotau1 = tau1;
            
            tau(r, d) = tau(r, d, index);
            tau1(r, d) = tempotau1(r, d, index);
           
            tau1_final(r, d) = tau1(r, d, index);
            Vchip_final(r, d) = Vchip(r, d, index);
            mu_final(r, d) = .1484 + .0061 .* Vchip_final(r, d);
            f_ang_final(r, d) = beta_n(r, d, index);
            %p0(r, d) = ((3.2 .* cosd(eta_s(r, d)) .* (cosd(f_ang_final(r, d))) .^ 2) ./ (cosd(eta_s(r, d)) .* sind(2 .* (s_ang_final(r, d) + f_ang_final(r, d) - alpha_n)))) .* s_stress_final(r, d);
            Lcf(r, d) = (chip(r, d) .* (zeta + 2) .* sind(s_ang_final(r, d) + f_ang_final(r, d) - alpha_n)) ./ (2 .* sind(s_ang_final(r, d)) .* cosd(f_ang_final(r, d)) .* cosd(eta_c));
            tau1Overp0f(r, d) = (Lcf(r, d) .* sind(s_ang_final(r, d)) .* cosd(eta_c) .* cosd(s_ang_final(r, d) + f_ang_final(r, d) - alpha_n)) ./ (chip(r, d) .* (zeta + 1) .* cosd(eta_s(r, d)) .* cosd(f_ang_final(r, d)));
            p00(r, d) = tau1_final(r, d) ./ tau1Overp0f(r, d);

            p0(r, d) = ((3.2 * cosd(eta_s(r, d)) * (cosd(f_ang_final(r, d))) ^ 2) / (cosd(eta_s(r, d)) * sind(2 * (s_ang_final(r, d) + f_ang_final(r, d) - alpha_n)))) * tau(r, d);
                        
            kt(r, d) = ktc{r, d}(index);
            kn(r, d) = kfc{r, d}(index);
            kr(r, d) = krc{r, d}(index);

            %% Ploughing force prerequisites
        
            % Critical depth of cut, hcr [mm]
            h_cr(r, d) = gs(r, d, 5) .* 1e-3 .* (1 - cosd(gs(r, d, 8)));
            
            % Total contact length on flank face
            l1(r, d) = deg2rad(gs(r, d, 8)) .* gs(r, d, 5) .* 1e-3;
            l2(r, d) = deg2rad(gs(r, d, 3)) .* gs(r, d, 5) .* 1e-3;
            l3(r, d) = ((cosd(gs(r, d, 3)) - cosd(gs(r, d, 8))) ./ sind(gs(r, d, 3))) .* gs(r, d, 5) .* 1e-3;
            lf(r, d) = l1(r, d) + l2(r, d) + l3(r, d);
            x_f{r, d} = linspace(0, lf(r, d), 100);
            a(r, d) = p0(r, d) ./ ((2 .* lf(r, d) .* l1(r, d)) - (lf(r, d) .^ 2));
            b(r, d) = - 2 .* a(r, d) .* l1(r, d);
            c(r, d) = p0(r, d);
            px{r, d} = (a(r, d) .* (x_f{r, d} .^ 2)) + (b(r, d) .* x_f{r, d}) + c(r, d);
            lst{r, d} = mu_final(r, d) .* px{r, d} - tau1(r, d);
        
            % This section is responsible for constructing the tau function at the third zone
            signer = @(vec) sign(vec);
            indexer = @(vec) find(diff(sign(vec)) ~= 0);
            signz{r, d} = signer(lst{r, d});
            interindi{r, d} = indexer(lst{r, d});
            
            if interindi{r, d}
                for ind = 1 : numel(interindi{r, d})
                    for int = interindi{r, d}(ind)
                        tau_f{r, d} = heaviside(-signz{r, d}) .* (mu_final(r, d) .* px{r, d}) + heaviside(signz{r, d}) .* (tau1(r, d));
                    end
                end
            else
                tau_f{r, d} = heaviside(-signz{r, d}) .* (mu_final(r, d) .* px{r, d}) + heaviside(signz{r, d}) .* (tau1(r, d));
            end
        end
    end

    time = ['The elapsed time for computing the force coefficients is ', num2str(toc(t)), ' seconds'];
end