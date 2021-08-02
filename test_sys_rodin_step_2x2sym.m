%% Test sys_rodin_step_2x2sym.m

sys_rodin_step_2x2sym

% Check models are identical
t_test = 1*(0:20)';
U_test = reshape(3*idinput(21*4), 21, 4);
U_test(:, 4) = 1;
[y_test, t_test] = lsim(Gpss, U_test, t_test);
[y_test2, t_test] = lsim(Gpd, U_test, t_test);

% figure(100)
% subplot(2,1,1)
% plot(t_test, y_test(:,1), t_test, y_test2(:,1));
% grid on
% legend('Gpss', 'Gpd','Interpreter','none')
% subplot(2,1,2)
% plot(t_test, y_test(:,2), t_test, y_test2(:,2));
% grid on
% legend('Gpss', 'Gpd','Interpreter','none')

assert(all(abs(y_test - y_test2) < 0.5e-2, [1 2]))
