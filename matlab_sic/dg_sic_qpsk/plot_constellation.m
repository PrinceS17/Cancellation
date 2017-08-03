function y = plot_constellation(x, y, y_clean)
% plot the constellation of x, y, y_clean

figure
subplot(2,2,1),plot(x,'.');
xlabel('In phase'); ylabel('Quadrature'); title('x');
title('x'); xlim([-2 2]); ylim([-2 2]);

subplot(2,2,2),plot(y,'.');
xlabel('In phase'); ylabel('Quadrature'); title('y');

subplot(2,2,3),plot(y_clean,'.');
xlabel('In phase'); ylabel('Quadrature');title('y_{clean}');

end