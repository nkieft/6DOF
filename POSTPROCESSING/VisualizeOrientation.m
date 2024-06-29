close all; clear all; clc

load("BRUHMOMENT.mat");
q_data = q_data(1:30:end, :);
q0 = q_data(:, 1);
q1 = q_data(:, 2);
q2 = q_data(:, 3);
q3 = q_data(:, 4);

for i = 1:length(q_data)
    orientationx(:, i) = BodyToInertial([1; 0; 0], q0(i), q1(i), q2(i), q3(i));
    orientationy(:, i) = BodyToInertial([0; 1; 0], q0(i), q1(i), q2(i), q3(i));
    orientationz(:, i) = BodyToInertial([0; 0; 1], q0(i), q1(i), q2(i), q3(i));
end
j = 1;

for i = 1:length(q_data)
    
    plot3([0, orientationx(3, i)], [0, orientationx(2, i)], [0, orientationx(1, i)], 'r');
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    hold on
    plot3([0, orientationy(3, i)], [0, orientationy(2, i)], [0, orientationy(1, i)], 'g');
    plot3([0, orientationz(3, i)], [0, orientationz(2, i)], [0, orientationz(1, i)], 'b');
    hold off
    disp(i/length(q_data));

end
