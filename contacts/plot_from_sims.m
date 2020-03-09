data=load("session3_sims5_taup0.0_1e4.mat");
for i=1:length(data.sims) %loop on the different sims for the same tau_p
    disp(i)
    plot(data.sims(i))
end