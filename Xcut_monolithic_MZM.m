clear all
import com.comsol.model.*
import com.comsol.model.util.*
model=mphload('simple_Xcut_monolithic_MZM_structure_COMSOL60.mph');

[Q1, unit] = mphint2(model,'abs(es.normE)*(abs(ewfd.Ex)^2)','surface');
[Q2, unit2] = mphint2(model,'abs(ewfd.Ex^2)','surface');

% ---
% mphint2 function returns MxP matrix. M indicates the number of inner
% solutions. Inner solutions are based on single frequency/time/parameter.
% Outer solutions stem from time-dependent simulations/parametric
% sweeps/eigenvalue studies(modal shapes)/ frequency domain based
% solutions.
% Inner solutions: Single frequency operation
% Outer Solutions: Eigenvalue dependent values (effective refractive indices for
% each solution
% ---

% Parameters for simulation
[G,upper_,def,d]=mphevaluate(model, 'gap');
[V]=mphevaluate(model, 'V0');
lambda=mphevaluate(model, 'lambda');
n_e=mphevaluate(model,'ne');
n_o=mphevaluate(model,'no');
n_eff_struct = mpheval(model, 'real(ewfd.neff)');%De las simulaciones
n_eff = n_eff_struct.d1(:,1);
r_33=30.9e-12; %De la bibliografía
L = 0.01; % 1cm

% v_pi = mpheval(model, '(lambda*G*ewfd.neff) ./ (2*n_e^4*r_33.* (G / V .* (Q1./Q2)))');

overlap = G / V .* (Q1./Q2);

v_pi = (lambda*G.*(n_eff.')) ./ (2*n_e^4*r_33.* (G / V .* (Q1./Q2)) );

v_pi_L_valid = v_pi / L;



