#need function to add sampling error to all points

#need function to handle data

data = read.csv("Data/tri_summary.csv")

reformatted_data = data.frame(time = unique(data$Time),
                              Be = data$Mean[which(data$Bacteria=="EryR")],
                              Bt = data$Mean[which(data$Bacteria=="TetR")],
                              Bet = data$Mean[which(data$Bacteria=="DRP")])

plot_output(reformatted_data)


#parameter values:
fixed_theta = c(tre = 0,      #ery gene transduction probability
                trt = 0,
                Nmax = 10^9)     #maximum bacteria population size

parms = c(Ge = 1.5,         #growth parameter for ery resistant
          Gt = 1.5,
          Get = 1.4)

#initial values:

yinit_mass = c(Be = min(reformatted_data$Be),    #resistant to ery
               Bt = min(reformatted_data$Bt),    #resistant to tet
               Bet = min(reformatted_data$Bet))    #resistant to ery and tet

trace = mcmc_run(run_model = run_mass, data = reformatted_data[1:10,],
                 yinit = yinit_mass, fixed_theta = fixed_theta,
                 theta = parms, theta_sd = c(0.5, 0.5, 0.5),
                 n_iter = 5000, likelihood = likelihood_func)


fitted_parms = c(fixed_theta, Ge = mean(trace[,"Ge"]), 
                 Gt = mean(trace[,"Gt"]), Get = mean(trace[,"Get"]))
fitted_parms

fitted_parms = c(fixed_theta, parms)

res_fit = run_mass(seq(0,24,1), yinit_mass, fitted_parms)

plot(log10(reformatted_data$Be), type = "p", ylim = c(1, max(log10(reformatted_data))), pch=19, col = "blue")
points(log10(reformatted_data$Bt), pch = 19, col = "green")
points(log10(reformatted_data$Bet), pch = 19, col = "red")
lines(log10(res_fit$Be), col = "blue")
lines(log10(res_fit$Bt), col = "green")
lines(log10(res_fit$Bet), col = "red")



plot(reformatted_data$Be, type = "p", ylim = c(1, max(reformatted_data)), pch=19, col = "blue")
points(reformatted_data$Bt, pch = 19, col = "green")
points(reformatted_data$Bet, pch = 19, col = "red")
lines(res_fit$Be, col = "blue")
lines(res_fit$Bt, col = "green")
lines(res_fit$Bet, col = "red")

