% MATLAB Code of the GPB2
% Adapted from code provided by Shunyi Zhao in private email
% See this publication:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

function [x,P,u,Out_x,Out_P] = GPB2_estimation(x,P,y,Model,u)

    TP = Model.TP;
    M = size(TP,1);
    N = size(Model.A{1}, 1);  % dimension of state
    A = Model.A;
    B = Model.B;
    C = Model.C;
    Q = Model.Q;
    R = Model.R;

    for j = 1:M

        for i = 1:M

            % M^2 Kalman filter predictions
            Pred_x = A{j} * x(:,i);
            Pred_P = A{j} * P(:,:,i) * A{j}' + Q{j};

            % Calculate the mixing probabilities
            S = C{j}  * Pred_P * C{j}'+ R{j};
            Mix_u(i,j) = TP(i,j) * u(i) * normpdf(y, C{j} * Pred_x, sqrt(S));

            % Update the M^2 KF estimates and error covariances
            K = Pred_P * C{j}' * S^-1;
            Updx(:,i) = Pred_x + K *(y - C{j} * Pred_x);
            UpdP(:,:,i) = Pred_P - K * C{j} * Pred_P;

        end

        % Normalize the mode probabilities
        Upd_u = Mix_u(:,j) / sum(Mix_u(:,j));

        % Mix the estimates
        Mix_x(:,j) = sum(Updx .* repmat(Upd_u', N, 1), 2);
        Mix_P(:,:,j) = zeros(N, N);
        for i = 1:M
            summP = Upd_u(i) * (UpdP(:,:,i) + ...
                (Updx(:,i) - Mix_x(:,j)) * (Updx(:,i) - Mix_x(:,j))');
            Mix_P(:,:,j) =  Mix_P(:,:,j) + summP;
        end

        Storu(j) = sum(Mix_u(:,j));

    end
    x = Mix_x;
    P = Mix_P;
    u = Storu / sum(Storu);
    % Above is same as u = sum(Mix_u, 1) ./ sum(Mix_u, [1 2])

    % Output
    Out_x = sum(x .* repmat(u, N, 1), 2);
    % Above same as Out_x = sum(x .* u, 2)
    Out_P = zeros(N, N);
    for i = 1:M
        summP =  u(i) * (P(:,:,i) + (x(:,i) - Out_x) * (x(:,i) - Out_x)');
        Out_P =  Out_P + summP;   
    end
    disp("GPB2_estimation complete")

end