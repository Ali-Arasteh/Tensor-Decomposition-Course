%% Question 1
img = im2double(imread('cameraman.tif'));
centered_img = img-mean(img);
normalized_img = normalize(img);
figure();
imshow(img);
title('original image');
figure();
imshow(centered_img);
title('centered image');
figure();
imshow(normalized_img);
title('normalized image');
[U, S, V] = svd(img);
[cen_U, cen_S, cen_V] = svd(centered_img);
[norm_U, norm_S, norm_V] = svd(normalized_img);
%% Part 1
% original image
figure();
i = 0;
for k=[2, 4, 8, 16, 32, 64, 128, 256]
    i = i+1;
    subplot(2, 4, i);
    approximated_img = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';
    imshow(approximated_img);
    title(['k = ', num2str(k), ', Error = ', num2str(norm(img-approximated_img, 'fro')/norm(img, 'fro'))]);
end
% centered image
figure();
i = 0;
for k=[2, 4, 8, 16, 32, 64, 128, 256]
    i = i+1;
    subplot(2, 4, i);
    approximated_img = cen_U(:, 1:k)*cen_S(1:k, 1:k)*cen_V(:, 1:k)'+mean(img);
    imshow(approximated_img);
    title(['k = ', num2str(k), ', Error = ', num2str(norm(img-approximated_img, 'fro')/norm(img, 'fro'))]);
end
% normalized image
figure();
i = 0;
for k=[2, 4, 8, 16, 32, 64, 128, 256]
    i = i+1;
    subplot(2, 4, i);
    approximate_norm_img = norm_U(:, 1:k)*norm_S(1:k, 1:k)*norm_V(:, 1:k)';
    imshow(approximate_norm_img);
    title(['k = ', num2str(k), ', Error = ', num2str(norm(normalized_img-approximate_norm_img, 'fro')/norm(normalized_img, 'fro'))]);
end
%% Part 2
% original image
error_2_img = zeros(1, 256);
error_fro_img = zeros(1, 256);
for k=1:256
    approximated_img = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';
    error_2_img(k) = norm(img-approximated_img)/norm(img);
    error_fro_img(k) = norm(img-approximated_img, 'fro')/norm(img, 'fro');
end
figure();
plot(1:256, error_2_img, 'r');
hold on;
grid on;
plot(1:256, error_fro_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('k');
ylabel('error');
title("Error between original image and approximated image");
xlim([1, 256]);
ylim([0, 1]);
% centered image
error_2_centered_img = zeros(1, 256);
error_fro_centered_img = zeros(1, 256);
for k=1:256
    approximated_img = cen_U(:, 1:k)*cen_S(1:k, 1:k)*cen_V(:, 1:k)'+mean(img);
    error_2_centered_img(k) = norm(img-approximated_img)/norm(img);
    error_fro_centered_img(k) = norm(img-approximated_img, 'fro')/norm(img, 'fro');
end
figure();
plot(1:256, error_2_centered_img, 'r');
hold on;
grid on;
plot(1:256, error_fro_centered_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('k');
ylabel('error');
title("Error between original image and approximated centered image");
xlim([1, 256]);
ylim([0, 1]);
% normalized image
error_2_normalized_img = zeros(1, 256);
error_fro_normalized_img = zeros(1, 256);
for k=1:256
    approximated_normalized_img = norm_U(:, 1:k)*norm_S(1:k, 1:k)*norm_V(:, 1:k)';
    error_2_normalized_img(k) = norm(normalized_img-approximated_normalized_img)/norm(normalized_img);
    error_fro_normalized_img(k) = norm(normalized_img-approximated_normalized_img, 'fro')/norm(normalized_img, 'fro');
end
figure();
plot(1:256, error_2_normalized_img, 'r');
hold on;
grid on;
plot(1:256, error_fro_normalized_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('k');
ylabel('error');
title("Error between normalized image and approximated normalized image");
xlim([1, 256]);
ylim([0, 1]);
%% Part 3
% original image
figure();
plot((256*2+1)/(256*256)*(1:256), error_2_img, 'r');
hold on;
grid on;
plot((256*2+1)/(256*256)*(1:256), error_fro_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('compression rate');
ylabel('error');
title("Error between original image and approximated image");
xlim([0, 2]);
ylim([0, 1]);
% centered image
figure();
plot((256*2+1+1)/(256*256)*(1:256), error_2_centered_img, 'r');
hold on;
grid on;
plot((256*2+1+1)/(256*256)*(1:256), error_fro_centered_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('compression rate');
ylabel('error');
title("Error between original image and approximated centered image");
xlim([0, 2]);
ylim([0, 1]);
% normalized image
figure();
plot((256*2+1)/(256*256)*(1:256), error_2_normalized_img, 'r');
hold on;
grid on;
plot((256*2+1)/(256*256)*(1:256), error_fro_normalized_img, 'g');
legend('|A - A_k|_2', '|A - A_k|_F');
xlabel('compression rate');
ylabel('error');
title("Error between normalized image and approximated normalized image");
xlim([0, 2]);
ylim([0, 1]);
%% Question 2
load('EEGdata.mat');
Xorg = Xorg';
Xnoise = Xnoise';
average_Xorg = mean(Xorg);
average_Xnoise = mean(Xnoise);
[Unoise, Snoise, Vnoise] = svd(Xnoise-average_Xnoise, 'econ');
% sources
source = Unoise*Snoise;
for k=1:32
    subplot(4, 8, k);
    plot(source(:, k));
    title(['Channel ', num2str(k)]);
end
% best Source
error_fro_one_signal = zeros(1, 32);
for k=1:32
    approximated_signal = Unoise(:, k)*Snoise(k, k)*Vnoise(:, k)'+average_Xnoise;
    error_fro_one_signal(k) = norm(Xorg-approximated_signal, 'fro');
end
figure();
stem(1:32, error_fro_one_signal);
grid on;
xlabel('k');
ylabel('|A - A_k one|_F');
title("Frobenius error between reconstructed image with source k and original image");
% k-best sources
error_fro_signal = zeros(1, 32);
for k=1:32
    approximated_signal = Unoise(:, 1:k)*Snoise(1:k, 1:k)*Vnoise(:, 1:k)'+average_Xnoise;
    error_fro_signal(k) = norm(Xorg-approximated_signal, 'fro');
end
figure();
stem(1:32, error_fro_signal);
grid on;
xlabel('k');
ylabel('|A - A_k|_F');
title("Frobenius error between reconstructed image with k-first sources and original image");
% best reconstructed signals
approximated_signal = Unoise(:, 1)*Snoise(1, 1)*Vnoise(:, 1)'+average_Xnoise;
best_error = norm(Xorg-approximated_signal, 'fro');
figure();
for k=1:32
    subplot(4, 8, k);
    plot(approximated_signal(:, k));
    title(['Channel ', num2str(k)]);
end
sgtitle(['Best recounstructed signals, Error = ', num2str(best_error)]);
figure();
for k=4:4:32
    subplot(2, 4, ceil(k/4));
    plot(approximated_signal(:, k), 'r');
    hold on;
    plot(Xorg(:, k), 'g');
    xlabel('t');
    ylabel('amplitude');
    legend('best denoised Signal', 'original Signal');
    title(['Channel ', num2str(k)]);
end
%% Question 3
%% Part 1
load('PCAdata.mat');
PCA_data = PCAdata';
average_data = mean(PCA_data);
figure();
subplot(3, 2, [1, 3, 5]);
scatter3(PCA_data(:, 1), PCA_data(:, 2), PCA_data(:, 3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('PCA data');
subplot(3, 2, 2);
scatter3(PCA_data(:, 1), PCA_data(:, 2), PCA_data(:, 3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('x-y plane');
view(0, 90);
subplot(3, 2, 4);
scatter3(PCA_data(:, 1), PCA_data(:, 2), PCA_data(:, 3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('x-z plane');
view(0, 0);
subplot(3, 2, 6);
scatter3(PCA_data(:, 1), PCA_data(:, 2), PCA_data(:, 3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('y-z plane');
view(90, 0);
figure();
subplot(3, 2, [1, 3, 5]);
scatter3(PCA_data(:, 1)-average_data(1), PCA_data(:, 2)-average_data(2), PCA_data(:, 3)-average_data(3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('centered PCA data');
subplot(3, 2, 2);
scatter3(PCA_data(:, 1)-average_data(1), PCA_data(:, 2)-average_data(2), PCA_data(:, 3)-average_data(3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('centered x-y plane');
view(0, 90);
subplot(3, 2, 4);
scatter3(PCA_data(:, 1)-average_data(1), PCA_data(:, 2)-average_data(2), PCA_data(:, 3)-average_data(3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('centered x-z plane');
view(0, 0);
subplot(3, 2, 6);
scatter3(PCA_data(:, 1)-average_data(1), PCA_data(:, 2)-average_data(2), PCA_data(:, 3)-average_data(3), 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g', 'LineWidth', 0.1);
xlabel('x');
ylabel('y');
zlabel('z');
title('centered y-z plane');
view(90, 0);
x_tautness = norm(PCA_data(:, 1));
y_tautness = norm(PCA_data(:, 2));
z_tautness = norm(PCA_data(:, 3));
centered_x_tautness = norm(PCA_data(:, 1)-average_data(1));
centered_y_tautness = norm(PCA_data(:, 2)-average_data(2));
centered_z_tautness = norm(PCA_data(:, 3)-average_data(3));
%% Part 2
[PCA_U, PCA_S, PCA_V] = svd(PCA_data-average_data, 'econ');
new_tautness_SVD = PCA_S;
new_variance_SVD = PCA_S.^2/1000;
new_direction_SVD = PCA_V;
new_data_SVD = PCA_U*PCA_S;
%% Part 3
[coeff, score, latent] = pca(PCA_data);
new_tautness_PCA = sqrt(latent*1000);
new_variance_PCA = latent;
new_direction_PCA = coeff;
new_data_PCA = score;
%% Question 4
load('house_dataset.mat');
house_inputs = houseInputs';
house_targets = houseTargets';
[house_U, house_S, house_V] = svd(normalize(house_inputs), 'econ');
source = house_U*house_S;
covariance = zeros(1, 13);
for i=1:13
    covariance(i) = (dot(house_targets, source(:, i))-mean(house_targets)*mean(source(:, i)))/(norm(house_targets)*norm(source(:, i)));
end
figure();
subplot(1, 2, 1);
scatter(source(:, 1), house_targets, 5, 'filled');
xlabel('first correlated source');
ylabel('house price');
subplot(1, 2, 2);
scatter(source(:, 3), house_targets, 5, 'filled');
xlabel('second correlated source');
ylabel('house price');
figure();
scatter(source(:, 1), source(:, 3), 10, house_targets, 'filled');
xlabel('first correlated source');
ylabel('second correlated source');
colorbar;