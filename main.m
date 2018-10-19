clear; clear all;
% set model variables for p: patients and c: controls
str_loop_p = [0.0 0.1 0.2]; str_loop_c = 0.;
str_noise = 0.2;

% record no. rejections for each value of str_loop_p
cnt_vec = zeros(1,length(str_loop_p));

parfor i = 1:length(str_loop_p)
    cnt = 0;
    for j=1:100 % no. experiments
        if mod(j,10) == 0
            disp(j);
        end
        res = experiment(str_loop_p(i), str_loop_c, str_noise);
        if res == 1
            cnt = cnt+1;
        end
    end
    cnt_vec(i) = cnt;
end

disp(cnt_vec)