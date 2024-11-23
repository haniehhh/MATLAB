clc;
clear;
% close all;
%%

material_inx = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% material = Sb2S3  -> material_inx = 1 %%%%
%%%% material = Sb2Se3 -> material_inx = 3 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch material_inx
    case 1
        material_name = "Sb_{2}S_{3}";
        filename = "1a";
    case 3
        material_name = "Sb_{2}Se_{3}";
        filename = "3b";
end

% if material_inx == 1
%      material_name = "Sb_{2}S_{3}";
%         filename = "1a";
% elseif material_inx == 3
%     material_name = "Sb_{2}Se_{3}";
%         filename = "3b";
% else
%     disp("wrong index")
% end



file = readtable(".\data for py\Diffuse\"+filename+" - R - NIR = 20nm - Diffuse.xlsx");


wavelength = file(:,1).Variables;
reflectance = file(:,2).Variables;

energy_eV = 1240./wavelength;
reflectance_normalized = reflectance/100;

energy_range = (energy_eV >1 ) & (energy_eV <3);

F = (1 - reflectance_normalized) .^ 2 ./ (2 * reflectance_normalized);
A = (energy_eV .* F) .^ 2; % to the power of 2, in the case of direct bandgap
% to the power of 0.5, in the case of indirect bandgap


% figure
% plot(energy_eV, F)
% xlim([1 2])
% title("$F(R_{\infty})$ = $\frac{ {(1-R_{\infty})}^2 }{ 2R_{\infty} }$", 'Interpreter', 'latex')
% xlabel("energy (eV)")
% ylabel("F Function")


% Find bandgap with smoothing the A function
[bandgap, sharpest_slope] = fit_findBG(A, energy_eV);

figure
plot(energy_eV, A, 'LineWidth',1.5)
xlim([1 2])
title("Kubelka-Munk Function")
xlabel("energy (eV)")
ylabel("{(F(R_{\infty})h{\nu})}^2")

hold on
lin_x = linspace(bandgap, bandgap+0.08, 20);
lin_y = sharpest_slope * (lin_x - bandgap);

plot(lin_x, lin_y, 'LineStyle','--', 'LineWidth', 2, 'Color', 'red' )

legend(material_name, "Fitted Line")

hold off



% bandgap_ev = py_findBG(A, energy_eV);

%%%%%%%%%%%%%%%%%%%%%%% Find Bandgap by fitting %%%%%%%%%%%%%%%%%%%%%%
function [bandgap, sharpest_slope] = fit_findBG(A, energy_eV)
    energy_range = (energy_eV >1 )& (energy_eV <2);

    energy_short = energy_eV(energy_range);
    A_short = A(energy_range);

    [f, gof] = fit(energy_short, A_short, 'gompertz');
    
    
    figure
    plot(energy_short, A_short, energy_short, f(energy_short))
    xlim([1 2])
    

    sharpest_grad = max( gradient( f(energy_short)));
    
    indx = find(gradient(f(energy_short))==sharpest_grad);

    sharpest_slope = sharpest_grad/(energy_short(indx) - energy_short(indx-1));
    x1 = energy_short(indx);
    y1 = A_short(indx);
    
    bandgap = - y1 / sharpest_slope + x1;


end



%%%%%%%%%%%%%%%%%%%%%%% Python code interpretation %%%%%%%%%%%%%%%%%%%%%%%
function [bandgap_ev] = py_findBG(A, energy_eV)

    energy_range = (energy_eV >1 )& (energy_eV <3);
    ds_factor = 1;
    A_ds = downsample(A,ds_factor);
    energy_ds = downsample(energy_eV, ds_factor);
    
    grad = gradient(A_ds);
    max_slope_indx = find(gradient(A) == max(gradient(A)) );
    max_grad = max(grad)/(energy_eV(max_slope_indx)-energy_eV(max_slope_indx-1));
    
    max_slope_energy = energy_eV(max_slope_indx);
    max_slope_A = A(max_slope_indx);
    
    bandgap_ev = - max_slope_A/max_grad + max_slope_energy;
    
    % figure
    % plot(energy_eV, gradient(A))
    
    energy_range_filtered = energy_eV(energy_range);
    A_range_filtered = A(energy_range);
    
    max_slope_index = find(gradient(A_range_filtered) == max(gradient(A_range_filtered)));
    window_size = 50; % Amount of points taken in the straight line
    max_r_squared = 0;
    best_window_start = 0;
    
    for i= max_slope_index : length(A_range_filtered) - window_size
        window_end = i + window_size;
        window_energy = energy_range_filtered(i:window_end);
        window_A = A_range_filtered(i:window_end);
    
        fitmit = fit(window_energy, window_A, 'poly1');
        slope = fitmit.p1;
        intercept = fitmit.p2;
        fitted = slope * window_energy + intercept;
    
        r_squared = corr(fitted, window_A);
    
        if r_squared > max_r_squared
            max_r_squared = r_squared;
            best_window_start = i;
        end
    end
    
    best_region_energy = energy_range_filtered(best_window_start:best_window_start + window_size);
    best_region_A = A_range_filtered(best_window_start:best_window_start + window_size);
    
    fitmit = fit(best_region_energy, best_region_A, 'poly1');
    slope = fitmit.p1;
    intercept = fitmit.p2;
    band_gap_energy = -intercept / slope ;
end


