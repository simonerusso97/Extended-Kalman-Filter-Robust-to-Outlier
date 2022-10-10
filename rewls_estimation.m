function res = rewls_estimation(z_pos, z_phi, N_gps, N_compass)
   
    %stima ls
    %x = fmincon(@(s)obj_fun(s, z_pos(1,:),[], N_gps ), 0,[],[],[],[],lb,ub);
    %y = fmincon(@(s)obj_fun(s, z_pos(2,:),[], N_gps ), 0,[],[],[],[],lb,ub);
    %phi = fmincon(@(s)obj_fun(s, z_phi, [], N_compass ), 0,[],[],[],[],lb,ub);
%     x = fminsearch(@(s)obj_fun(s, z_pos(1,:),[], N_gps ),0);
%     y = fminsearch(@(s)obj_fun(s, z_pos(2,:),[], N_gps ),0);
%     phi = fminsearch(@(s)obj_fun(s, z_phi,[], N_gps ),0);
    

    for i=1:length(z_pos(1,:))
        H_x(i,1) = 1;
    end
    H_y = H_x;

    x = (H_x' * H_x)^-1 * H_x' *  z_pos(1,:)';
    y = (H_y' * H_y)^-1 * H_y' *  z_pos(2,:)';


    for i=1:length(z_phi(1,:))
        H_phi(i,1) = 1;
    end

    phi = (H_phi' * H_phi)^-1 * H_phi' + z_phi';




    %calcolo residui e sigma robusta agli outlier
    c = 1.4826;

    r_pos(1,:) = z_pos(1,:) - x;
    r_pos(2,:) = z_pos(2,:) - y;

    r_phi = z_phi(:) - phi;

    r_pos(1,:) = sort(r_pos(1,:));
    r_pos(2,:) = sort(r_pos(2,:));
    r_phi(:) = sort(r_phi(:));

    mu(1) = median(r_pos(1,:));
    mu(2) = median(r_pos(2,:));
    mu(3) = median(r_phi(:));

    madd(1) = median(abs(r_pos(1,:)-mu(1)));
    madd(2) = median(abs(r_pos(2,:)-mu(2)));
    madd(3) = median(abs(r_phi(:)-mu(3)));
    

    %mad(X,1) computes Y based on medians, i.e. MEDIAN(ABS(X-MEDIAN(X))
    sigma_x_hat = c*madd(1);%c * mad(r_pos(1,:),1);
    sigma_y_hat = c*madd(2);%c * mad(r_pos(2,:),1);
    sigma_phi_hat = c*madd(1);%c * mad(r_phi,1);

    %calcolo i pesi w
    
    w_x_star = compute_w(r_pos(1,:), sigma_x_hat);

    w_y_star = compute_w(r_pos(2,:), sigma_y_hat);

    w_phi_star = compute_w(r_phi, sigma_phi_hat);

    %calcolo sigma_star

    num = 0;
    den = -1;
    for i = 1:N_gps
        num = num + w_x_star(i)*r_pos(1,i)^2;
    den = den + w_x_star(i);
    end

    sigma_x_star =  sqrt ( num/den );

    num = 0;
    den = -1;
    for i = 1:N_gps
        num = num + w_y_star(i)*r_pos(2,i)^2;
    den = den + w_y_star(i);
    end
    

    sigma_y_star =  sqrt ( num/den );

    num = 0;
    den = -1;
    for i = 1:N_compass
        num = num + w_phi_star(i)*r_phi(i)^2;
    den = den + w_phi_star(i);
    end

    sigma_phi_star =  sqrt ( num/den );

    %calcolo w_sword
    w_x_sword = compute_w(r_pos(1,:), sigma_x_star);
    w_y_sword = compute_w(r_pos(2,:), sigma_y_star);
    w_phi_sword = compute_w(r_phi, sigma_phi_star);

    %stima RLS

    R_x = diag(w_x_sword);
    R_y = diag(w_y_sword);
    R_phi = diag(w_phi_sword);

    x_rls = (H_x' * R_x^-1 * H_x)^-1 * H_x' * R_x^-1 * z_pos(1,:)'; %fmincon(@(s)obj_fun(s, z_pos(1,:), w_x_sword, N_gps ), 0,[],[],[],[],lb,ub);
    y_rls = (H_y' * R_y^-1 * H_y)^-1 * H_y' * R_y^-1 * z_pos(2,:)'; %fmincon(@(s)obj_fun(s, z_pos(2,:), w_y_sword, N_gps ), 0,[],[],[],[],lb,ub);
    phi_rls =(H_phi' * R_phi^-1 * H_phi)^-1 * H_phi' * R_phi^-1 * z_phi'; %fmincon(@(s)obj_fun(s, z_phi, w_phi_sword, N_compass ), 0,[],[],[],[],lb,ub);

    res(1,:) = [x_rls y_rls phi_rls];
    res(2, :) = [sigma_x_star sigma_y_star sigma_phi_star];


end