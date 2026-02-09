// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;


// Gamma_func2 = function(e, x, d) {
//   n = nrow(x)
//   p = ncol(x)
//   result = matrix(0,n,p)
//   emat = outer(e,e,"<") # indicator function
//   term1 = matrix(0,n,p)
//   for (i in 1:n) {
//     term1[i,] = colSums(emat[i,]*d[i]*(matrix(x[i,],n,p,byrow=T)-x))
//   }
//   term2 = matrix(0,n,p)
//     for (i in 1:n) {
//       xi = x[i,]
//       delta_i = d[i]
//       for (j in 1:n) {
//         if (e[i] >= e[j]) {
//           eej = which(e>=e[j])
//           numer = colSums(matrix(xi,length(eej),p,byrow=T)-x[eej,drop=FALSE])
//           denom = length(eej)
//           if (denom > 0) {
//             term2[i,] = term2[i,]+(d[j])*(numer/denom)
//           }
//         }
//       }
//     }
//     result = term1-term2
//     cov(result)
// }

// [[Rcpp::export]]
arma::mat Gamma_func2_cpp(const arma::vec& e, const arma::mat& x, const arma::vec& d) {
  int n = x.n_rows;
  int p = x.n_cols;

  arma::mat term1(n,p,fill::zeros);
  arma::mat term2(n,p,fill::zeros);

  arma::umat emat = arma::zeros<arma::umat>(n,n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      emat(i,j) = (e[i] < e[j]);
    }
  }
  for (int i = 0; i < n; i++){
    arma::rowvec xi = x.row(i);
    arma::mat diff = arma::repmat(xi, n, 1) - x; // matrix(x[i,],n,p,byrow=T)-x
    arma::rowvec sumvec = arma::sum(arma::repmat(emat.row(i), p, 1).t() % diff, 0); // emat[i,]* (matrix)
    term1.row(i) = d[i] * sumvec; //d[i]*~
  }

  for (int i = 0; i < n; i++) {
    arma::rowvec xi = x.row(i);
    for (int j = 0; j < n; j++){
      if (e[i] >= e[j]) {
        arma::uvec eej = arma::find(e >= e[j]); // eej = which(e >= e[j])
        int denom = eej.n_elem;
        if (denom > 0) {
          arma::mat xeej = x.rows(eej); // x[eej, drop=FALSE]
          arma::mat diff = arma::repmat(xi, denom, 1) - xeej; // matrix(xi,length(eej),p,byrow=T)
          arma::rowvec numer = arma::sum(diff, 0); // colSums()
          term2.row(i) += d[j] * (numer / denom); //  term2[i,] = term2[i,]+(d[j])*(numer/denom)
        }
      }
    }
  }
  arma::mat result = term1-term2;
  return arma::cov(result);
}


// un_est_func = function(ei,ej,xij,di,si,z=NULL) {
//   n = sqrt(length(si))
//   ind_i = rep(1:n,each=n)
//   if (is.null(z)) z = rep(1,n)
//     zi = z[ind_i]
//   Phi = zi*(ei<ej) # indicated function
//
//   drop(t(di*xij) %*% Phi)/n
// }

// [[Rcpp::export]]
arma::vec un_est_func_cpp(const arma::vec& ei,
                          const arma::vec& ej,
                          const arma::mat& xij,
                          const arma::vec& di,
                          const arma::vec& z) {
  int n = std::sqrt(ei.n_elem);
  arma::vec zi(ei.n_elem);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      zi[i*n + j] = z[i];
    }
  }

  arma::vec Phi = zi % (ei < ej);
  arma::vec tmp = di % Phi;
  arma::vec out = (xij.t() * tmp) / n;

  return out;
}



// est_func = function(ei,ej,xij,di,si,z=NULL) {
//   n = sqrt(length(si))
//   ind_i = rep(1:n,each=n)
//   if (is.null(z)) z = rep(1,n)
//     zi = z[ind_i]
//   Phi = zi*drop(pnorm((ej-ei)/si)) # standard normal distribution function
//
//   drop(t(di*xij) %*% Phi)/n
// } # first derivative func

// [[Rcpp::export]]
arma::vec est_func_cpp(const arma::vec& ei,
                   const arma::vec& ej,
                   const arma::mat& xij,
                   const arma::vec& di,
                   const arma::vec& si) {
  int n = std::sqrt(ei.n_elem);
  arma::vec zi(ei.n_elem);
  arma::vec z = arma::ones<arma::vec>(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      zi[i * n + j] = z[i];
    }
  }
  arma::vec ratio = (ej - ei) / si;
  arma::vec Phi(ratio.n_elem);

  for (uword k = 0; k < ratio.n_elem; ++k) {
    Phi[k] = R::pnorm5(ratio[k], 0.0, 1.0, 1, 0 );
    }
  arma::vec weight = zi % Phi;
  arma::vec out = (xij.t() * (di % weight)) / n;

  return out;
}

// Gamma_func = function(ei,ej,xij,di,si,B=200) {
//   n = sqrt(length(si))
//   Shat = t(replicate(B,{
//     un_est_func(ei,ej,xij,di,si,z=rexp(n))*sqrt(n)
//   }))
//
//   cov(Shat)
// }

// [[Rcpp::export]]
arma::mat Gamma_func_cpp(const arma::vec& ei,
                         const arma::vec& ej,
                         const arma::mat& xij,
                         const arma::vec& di,
                         const arma::vec& si,
                         const int B = 200) {
  int n = std::sqrt(si.n_elem);
  int p = xij.n_cols;
  arma::mat Shat(B, p);

  for (int b = 0; b < B; b++) {
    arma::vec z(n);
    for (int i = 0 ; i < n; i++) {
      z[i] = R::rexp(1.0);
    }
    arma::vec score = un_est_func_cpp(ei, ej, xij, di, si);
    Shat.row(b) = (score.t() * std::sqrt((double)n));
  }
  arma::rowvec mean_row = arma::mean(Shat, 0);
  arma::mat centered = Shat.each_row() - mean_row;
  arma::mat Gamma = (centered.t() * centered) / (B-1) ;

  return Gamma;
}

// up_Amat = function(ei,ej,xij,di,si) {
//   phi = drop(dnorm((ej-ei)/si))
//
//   crossprod(di*xij*phi/si,xij)/n
// } # hessian

// [[Rcpp::export]]
arma::mat up_Amat_cpp(const arma::vec& ei,
                  const arma::vec& ej,
                  const arma::mat& xij,
                  const arma::vec& di,
                  const arma::vec& si) {
  int n = std::sqrt(ei.n_elem);
  arma::vec ratio = (ej - ei) / si;
  arma::vec phi(ratio.n_elem);

  for (uword k = 0; k < ratio.n_elem; ++k) {
    phi[k] = R::dnorm4(ratio[k], 0.0, 1.0, 0);
  }
  arma::vec weight = di % phi / si;
  arma::mat A = (xij.each_col() % weight).t() *  xij/n;

  return A;
}

// Sigma_func = function(ei,ej,xij,di,Gamma,si) {
//   A = up_Amat(ei,ej,xij,di,si)
//   Ainv = solve(A)
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_func_cpp(const arma::vec& ei,
                     const arma::vec& ej,
                     const arma::mat& xij,
                     const arma::vec& di,
                     const arma::mat& Gamma,
                     const arma::vec& si) {
  int n = std::sqrt(ei.n_elem);
  arma::mat A = up_Amat_cpp(ei, ej, xij, di, si);
  arma::mat Ainv;
  bool ok = arma::inv_sympd(Ainv, A);
  if(!ok) {
    Ainv = arma::inv(A);
  }
  arma::mat Sigma = Ainv * Gamma * Ainv / n;
  return Sigma;
}

// Sigma_func_pen = function(ei,ej,xij,di,Gamma,si,tilde_beta, lambda) {
//   A = up_Amat(ei,ej,xij,di,si) + diag((2*lambda) / (tilde_beta^2 + 1e-05))
//   Ainv = solve(A)
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_func_pen_cpp(const arma::vec& ei,
                         const arma::vec& ej,
                         const arma::mat& xij,
                         const arma::vec& di,
                         const arma::mat& Gamma,
                         const arma::vec& si,
                         const arma::vec& tilde_beta,
                         double lambda) {
  int n = std::sqrt(ei.n_elem);
  int p = tilde_beta.n_elem;

  arma::mat A = up_Amat_cpp(ei, ej, xij, di, si);
  arma::vec pen = (2 * lambda) / (arma::square(tilde_beta) + 1e-05);
  A.diag() += pen;
  arma::mat Ainv;
  bool ok = arma::inv_sympd(Ainv, A);
  if (!ok) {
    Ainv = arma::inv(A);
  }
  arma::mat Sigma = Ainv * Gamma * Ainv / n;
  return Sigma;
}

// sifunc = function(xij, Sigma) {
//   si = pmax(0.0001, rowSums((xij %*% Sigma) * xij))
//
//   sqrt(si)
// } # sifunc = ((Xi-Xj)*sigma*(Xi-Xj))^(1/ 2)

// [[Rcpp::export]]
arma::vec sifunc_cpp(const arma::mat& xij,
                 const arma::mat& Sigma) {
  int npair = xij.n_rows;
  int p = xij.n_cols;

  arma::vec si(npair);
  arma::mat M = xij * Sigma;
  for (int i = 0; i < npair; i++) {
    double v = arma::dot(M.row(i), xij.row(i));
    si[i] = std::max(0.0001, v);
  }
  return arma::sqrt(si);
}

// is_aft = function(data,se="CF") {
//   y = dt$y; d = dt$d; x = as.matrix(dt[,-(1:2)])
//     n = nrow(x)
//     p = ncol(data[,-c(1,2)])
//     Sigma = diag(p)/n
//   old_beta = lm(log(y) ~ x, weights = d)$coef[-1]
//   ind_i = rep(1:n,each=n)
//     ind_j = rep(1:n,times=n)
//     xi = x[ind_i,]
//   xj = x[ind_j,]
//   xij = xi-xj
//   yi = y[ind_i]
//   yj = y[ind_j]
//   di = d[ind_i]
//   ei = log(yi) - xi%*%old_beta
//   ej = log(yj) - xj%*%old_beta
//   si = sifunc(xij,Sigma)
//
//     maxit = 100; iter = 0; err = 10; tol = 1e-2
//     while (iter < maxit & tol < err) {
//       A = up_Amat(ei,ej,xij,di,si)
//       S = est_func(ei,ej,xij,di,si)
//       new_beta = drop(old_beta - S %*% solve(A)) # newton raphson method
//       err = norm(new_beta - old_beta,"2") # check for convergence of beta
//       old_beta = new_beta
//       e = log(y) - drop(x%*%old_beta)
//         ei = log(yi) - xi%*%old_beta
//       ej = log(yj) - xj%*%old_beta
//       if (se == "ZL") Gamma = Gamma_func(ei,ej,xij,di,si)
//         if (se == "CF") Gamma = Gamma_func2(e,x,d)
//           Sigma = Sigma_func(ei,ej,xij,di,Gamma,si)
//           si = sifunc(xij,Sigma)
//           iter = iter + 1
//     }
//     obj = list(beta = new_beta,
//                hess = A,
//                grad = S)
// }

// [[Rcpp::export]]
arma::vec weighted_lm_cpp(const arma::vec& y,
                      const arma::mat& x,
                      const arma::vec& d) {
  arma::vec logy = arma::log(y);
  arma::mat W = arma::diagmat(d);
  arma::mat XtWX = x.t() * W * x;
  arma::vec XtWy = x.t() * W * logy;

  arma::vec beta = arma::solve(XtWX, XtWy);
  return beta;
}


// [[Rcpp::export]]
Rcpp::List is_aft_cpp(const arma::vec& y,
                  const arma::vec& d,
                  const arma::mat& x,
                  std::string se_type = "") {
  int n = x.n_rows;
  int p = x.n_cols;

  arma::mat Sigma = arma::eye(p, p) / n;
  arma::vec old_beta = weighted_lm_cpp(y, x, d);
  arma::uvec ind_i(n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ind_i[i*n + j] = i;
    }
  }
  arma::uvec ind_j(n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ind_j[i * n + j] = j;
    }
  }
  arma::mat xi = x.rows(ind_i);
  arma::mat xj = x.rows(ind_j);
  arma::mat xij = xi-xj;

  arma::vec yi = y.elem(ind_i);
  arma::vec yj = y.elem(ind_j);

  arma::vec di = d.elem(ind_i);
  arma::vec ei = arma::log(yi) - (xi * old_beta);
  arma::vec ej = arma::log(yj) - (xj * old_beta);

  arma::vec si = sifunc_cpp(xij, Sigma);
  int maxit = 100;
  double tol = 1e-02;
  double err = 10;
  int iter = 0;

  arma::mat A;
  arma::vec S;
  arma::vec new_beta;
  arma::mat Gamma;

  while(iter < maxit && err > tol) {
    A = up_Amat_cpp(ei, ej, xij, di, si);
    S = est_func_cpp(ei, ej, xij, di, si);
    new_beta = old_beta - arma::solve(A, S);

    err = arma::norm(new_beta - old_beta, 2);
    old_beta = new_beta;
    arma::vec e = arma::log(y) - (x * old_beta);
    ei = arma::log(yi) - (xi * old_beta);
    ej = arma::log(yj) - (xj * old_beta);

    if (se_type == "ZL") {
      Gamma = Gamma_func_cpp(ei, ej, xij, di, si);
    } else {
      Gamma = Gamma_func2_cpp(e, x, d);
    }
    Sigma = Sigma_func_cpp(ei, ej, xij, di, Gamma, si);
    si = sifunc_cpp(xij, Sigma);

    iter++;
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = new_beta,
    Rcpp::Named("hess") = A,
    Rcpp::Named("grad") = S
  );
}

// Sigma_func = function(ei,ej,xij,di,Gamma,si, stable_numb1) {
//   A = up_Amat(ei,ej,xij,di,si) + diag(stable_numb1, p)
//   Ainv = solve(A)
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_func_stable_cpp(const arma::vec& ei,
                                const arma::vec& ej,
                                const arma::mat& xij,
                                const arma::vec& di,
                                const arma::mat& Gamma,
                                const arma::vec& si,
                                double stable_numb1) {
  int p = xij.n_cols;
  arma::mat A = up_Amat_cpp(ei, ej, xij, di, si);
  A.diag() += stable_numb1;
  arma::mat Ainv = arma::inv(A);
  int n = std::sqrt((double)ei.n_elem);
  return Ainv * Gamma * Ainv / n;
}


// Sigma_func_pen = function(ei,ej,xij,di,Gamma,si,tilde_beta, lambda, stable_numb1) {
//   A = up_Amat(ei,ej,xij,di,si) + diag((2*lambda) / (tilde_beta^2 + 1e-05))
//   A = A + diag(stable_numb1, p)
//   Ainv = solve(A)
//   Ainv %*% Gamma %*% Ainv/n
// }

arma::mat Sigma_func_pen_stable_cpp(const arma::vec& ei,
                                    const arma::vec& ej,
                                    const arma::mat& xij,
                                    const arma::vec& di,
                                    const arma::mat& Gamma,
                                    const arma::vec& si,
                                    const arma::vec& tilde_beta,
                                    double lambda,
                                    double stable_numb1) {
  int p = xij.n_cols;
  arma::mat A = up_Amat_cpp(ei, ej, xij, di, si);
  arma::vec pen = (2.0*lambda) / (arma::square(tilde_beta) + 1e-5);
  A.diag() += pen;
  A.diag() += stable_numb1;
  arma::mat Ainv = arma::inv(A);
  int n = std::sqrt((double)ei.n_elem);
  return Ainv * Gamma * Ainv / n;
}


// [[Rcpp::export]]
Rcpp::List is_aft_stable_cpp(const arma::vec& y,
                      const arma::vec& d,
                      const arma::mat& x,
                      const double stable_numb1,
                      std::string se_type = "") {
  int n = x.n_rows;
  int p = x.n_cols;

  arma::mat Sigma = arma::eye(p, p) / n;
  arma::vec old_beta = weighted_lm_cpp(y, x, d);
  arma::uvec ind_i(n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ind_i[i*n + j] = i;
    }
  }
  arma::uvec ind_j(n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ind_j[i * n + j] = j;
    }
  }
  arma::mat xi = x.rows(ind_i);
  arma::mat xj = x.rows(ind_j);
  arma::mat xij = xi-xj;

  arma::vec yi = y.elem(ind_i);
  arma::vec yj = y.elem(ind_j);

  arma::vec di = d.elem(ind_i);
  arma::vec ei = arma::log(yi) - (xi * old_beta);
  arma::vec ej = arma::log(yj) - (xj * old_beta);

  arma::vec si = sifunc_cpp(xij, Sigma);
  int maxit = 100;
  double tol = 1e-02;
  double err = 10;
  int iter = 0;

  arma::mat A;
  arma::vec S;
  arma::vec new_beta;
  arma::mat Gamma;

  while(iter < maxit && err > tol) {
    A = up_Amat_cpp(ei, ej, xij, di, si);
    S = est_func_cpp(ei, ej, xij, di, si);
    new_beta = old_beta - arma::solve(A, S);

    err = arma::norm(new_beta - old_beta, 2);
    old_beta = new_beta;
    arma::vec e = arma::log(y) - (x * old_beta);
    ei = arma::log(yi) - (xi * old_beta);
    ej = arma::log(yj) - (xj * old_beta);

    if (se_type == "ZL") {
      Gamma = Gamma_func_cpp(ei, ej, xij, di, si);
    } else {
      Gamma = Gamma_func2_cpp(e, x, d);
    }
    Sigma = Sigma_func_stable_cpp(ei, ej, xij, di, Gamma, si,stable_numb1);
    si = sifunc_cpp(xij, Sigma);

    iter++;
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = new_beta,
    Rcpp::Named("hess") = A,
    Rcpp::Named("grad") = S
  );
}
