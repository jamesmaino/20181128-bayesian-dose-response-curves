data {

  int <lower = 0> N; 
  int <lower = 0, upper = 1> y[N] ;
  int <lower = 1> pop[N];
  vector <lower = 0> [N] dose; 

}

parameters {

  vector [2] a;
  vector [2] b;
  
}

model {
     // Prior models for the unobserved parameters
    a ~ normal(0, 10);
    b ~ normal(1, 10);

    y ~ bernoulli_logit(a[pop] + (b[pop] .* log(dose)));

}
