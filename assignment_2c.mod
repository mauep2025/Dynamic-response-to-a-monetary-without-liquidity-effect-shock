// Basic New Keynesian Model (Gali) 2003

//-----------------------------------------------------------------//
//       1.  VARIABLES & PARAMETERS OF THE MODEL                   //
//-----------------------------------------------------------------//


var infl, y, ir, rn, a, ye, r, yn, dm, m, n;
varexo e_m e_a;

parameters mu, beta, sigma, rhom, phi, phina, lambda, kappa, theta, alpha, epsilon, rho, eta, gamma, nu; 

//--------------------------------------------------------------------//
//        2.PARAMETERS OF THE MODEL(CALIBRATION)                      //                
//--------------------------------------------------------------------//

alpha = 1/3;
phi = 1;
epsilon = 11;
beta = 0.99;
theta = .75;
sigma = 1;
rho = -log(beta);
rhom = .5;
eta = 4;
nu=1;
// define parameters that are functions of other parameters:
lambda = ((1-theta)*(1-beta*theta)/theta);
kappa = lambda*(sigma + phi);
phina =  (phi+1)/(sigma+phi);
mu = log(epsilon/(epsilon-1));
gamma = (nu-mu)/(sigma+phi);

//----------------------------------------------------------------//
//                3.EQUATIONS OF THE MODEL                        //
//----------------------------------------------------------------//

model;

/////////// Basic equations 

// New Keynesian Phillips Curve
infl = beta*infl(+1)+kappa*y;
// Dynamic IS
y = y(+1)-(1/sigma)*(ir-infl(+1)-rn);

// (Fisher or real interest rate)
r = ir-infl(+1);
// production function
ye = a+(1-alpha)*n;

////////// Auxiliary equations: 
// Natural Interest Rate
rn = rho+sigma*phina*(a(+1)-a);
// adhoc money (growth) demand
yn-eta*ir=m-y;
m(-1)=m+infl-dm;
// natural product and definition of output gap
yn = phina*a+gamma;
y = ye+yn;


// exogenous shocks
// monetary shock
dm = rhom*dm(-1)+e_m;

// non active productivity shock
a = 0*a(-1)+e_a;

end;

//-----------------------------------------------------------------------//
//                4.INITIAL VALUES FOR THE STEADY-STATE                  //
//-----------------------------------------------------------------------//


initval;
ir=0;

m = 0;
yn = 0;
ye = 0;
infl = 0;
a = 0;
r = 0;
y = 0;
ir = 0;
rn = 0;
dm = 0;
n = 0;
e_m = 0;
e_a = 0;

end;

steady;
check;

//-----------------------------------------------------------------------//
//             5.ST.DEV OF STOCHASTIC SHOCKS                             //
//-----------------------------------------------------------------------//

shocks;

var e_m; 
stderr 1;

var e_a; 
stderr 1;

end;

//-----------------------------------------------------------------------//
//            7.SIMULATION-RESULTS OF THE MODEL AND GRAPHS
//-----------------------------------------------------------------------//

stoch_simul( irf=25, order=1, nodisplay);

t=0:1:24;
figure(1);
title('Dynamic Responses to a monetary Shock');

subplot(4,1,1);
plot(t,4*infl_e_m,'-ob');
title('Inflation');

subplot(4,1,2);
plot(t,ye_e_m,'-ob');
title('Output');

subplot(4,1,3);
plot(t,4*r_e_m,'-ob');
title('Real Rate');

subplot(4,1,4);
plot(t,4*ir_e_m,'-ob');
title('Nominal Rate');
%axis([0 12 -10 10]);




