clear all;
close all;
clc;

rng(1);

%parametri automobile
l_r = 2.30; %m
l_f = 2.00; %m

%paramatri vari
Tc = 0.01; %passo di campionamento
T = 100; %running time
N = round(T/Tc); % n. sample
b = 0.15; %soglia outlier
N_gps = 1; %misurazioni al secondo del GPS
N_compass = 100; %misurazioni al secondo del compass sistem


%incertezze su s
sigma_x = 0.01; %m
sigma_y = 0.01; %m
sigma_phi = 0.1; %rad
Q = diag([sigma_x^2 sigma_y^2 sigma_phi^2]);


%incertezza sulla misura
sigma_z1 = 0.05; %m
sigma_z2 = 0.05; %m
sigma_z3= 0.0017; %rad
R = diag([sigma_z1^2 sigma_z2^2 sigma_z3^2]);

%vettore di stato iniziale reale
s_real(:,1) = [0 0 0]; %[m m rad]

%stima posozione iniziale

%acquisizione delle misure per 10 secondi a veicolo fermo
for t = 1:(10*N_gps)
    z_x(t) = s_real(1,1) + sigma_x *randn;
    z_y(t) = s_real(2,1) + sigma_y *randn;
end
for t = 1:(10*N_compass)
     z_phi(t) = s_real(3,1) + sigma_phi *randn;
end

%genero gli outlier
z_x(3) = z_x(3)+ 2;
z_y(1) = z_y(1)+ 10;
z_phi(34) = z_phi(34) + 25*randn;
z_phi(45) = z_phi(45) + 25*randn;


z_pos(1,:) = z_x;
z_pos(2,:) = z_y;

%qualcosa non va, si verifica un problema di tipo numerico
res = rewls_estimation(z_pos, z_phi, 10*N_gps, 10*N_compass);
res = [ 0.01, 0.02 0.015; 0.02 0.02 0.035];
%EKF 

%ingressi di controllo
u(1,1:N) = 10; %v in m/s;
u(2,1:N) = 0.0873; %delta in radianti

%inizzializzazione
s_ekf = res(1,:)';
P_ekf = diag(res(2,:).^2);

for k = 1:N
    %modello reale s
    s_real(:,k+1) = s_evolution_model(s_real(:,k), u(:,k), Tc, l_r, l_f)' + sqrt(Q) * randn(3,1);

    %predizione di misura di s

    z_pred(:,k+1)  =  predict_z(s_real(:,k+1), k, Tc, N_gps, N_compass);

    z(:, k+1) = z_pred(:,k+1) + sqrt(R) * randn(3,1);
    
    %generazione outlier
    contamination_x = 100*rand;
    contamination_y = 100*rand;
    contamination_phi = 100*rand;
    
    if mod( (k-1) , (N_gps/Tc) ) == 0
        if  contamination_x > 45 && contamination_x < 55
            z(1, k+1) = z(1,k+1) + 25*randn;
        end
        if  contamination_y > 45 && contamination_y < 55
            z(2, k+1) = z(2,k+1) + 25*randn;
        end
    end

    if  contamination_phi > 45 && contamination_phi < 55
            z(3, k+1) = z(3,k+1) + 25*randn;
    end


     %EKF1 PREDIZIONE DELLO STATO
    beta = atan( l_r * tan(u(2,k)) / ( l_r+l_f) );
    F = compute_F(Tc, u(1,k), s_ekf(3,k), beta);
    L = eye(3);

    P_pred = F * P_ekf * F' + L * Q * L';
    s_pred(:,k+1) = s_evolution_model(s_ekf(:,k), u(:,k), Tc, l_r, l_f)';

    %EKF1 CORREZIONE DELLO STATO

    clear ni;
    if not(isnan(z(1,k+1)))

        H = eye(3);
        M = eye(3);
        
        W = P_pred * H' * ( H * P_pred * H' + M * R * M' )^-1;
  
        ni(:) = z(:,k+1) - s_pred(:, k+1);
        
    else
        H = [0 0 1];
            
        M = [0 0 1];

        W = P_pred * H' * ( H * P_pred * H' + M * R * M' )^-1;

        ni = z(3,k+1) - s_pred(3, k+1);

    end

        P_ekf = ( eye(3) - W * H ) * P_pred;
        s_ekf(:,k+1) = s_pred(:,k+1) + W * ni(:) * min(1, b/ ( norm( W * ni(:) ) ) );

end

figure;
hold on;
plot(s_real(1,:), s_real(2,:),'Color', 'g', 'LineWidth', 1.3);
plot(s_ekf(1,:), s_ekf(2,:), 'r');
legend('traiettoria reale', 'traiettoria stimata EFK');

figure;
hold on;
plot(wrapTo360( rad2deg( s_real(3,:) ) ), 'Color', 'g', 'LineWidth', 1.3);
plot(wrapTo360( rad2deg( s_ekf(3,:) ) ), 'r');
legend('phi reale', 'phi stimato EKF');

figure;
plot( sqrt((s_real(1,:) - s_ekf(1,:) ).^2 + (s_real(2,:) - s_ekf(2,:) ).^2), 'b' );
legend("square error traiettoria");

figure;
plot( wrapTo360( rad2deg( sqrt( ( s_ekf(3,:) - s_real(3,:) ).^2 ) ) ) , 'b' );
legend("square error orientazione");

figure;
hold on
plot( wrapTo360( rad2deg( z(3,:) ) ) , 'b' );
plot( wrapTo360( rad2deg( s_real(3,:) ) ) , 'Color', 'g', 'LineWidth', 1.3 );
legend("orientation measure", "orientation real");

k=2:100:N;

figure;
hold on;
plot(k, s_real(1,k),'Color', 'g', 'LineWidth', 1.3);
plot(k, z(1,k), 'r');
legend('real x position ', 'x position measure');

figure;
hold on;
plot(k, s_real(2,k),'Color', 'g', 'LineWidth', 1.3);
plot(k, z(2,k), 'r');
legend('real y position ', 'y position measure');