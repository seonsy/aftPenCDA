// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;


// Gamma_func2p = function(dt, beta) {
//   n = length(unique(dt$id))
//   id = dt$id
//   m = as.numeric(table(id))
//   N = sum(m)
//   x = as.matrix(dt[, -(1:4)])
//   p = ncol(x)
//   y = pmax(dt$L, 1e-8)
//   d = dt$delta
//   v = ifelse(d == 1, log(y), dt$R) - x %*% beta
//   u = ifelse(d == 1, log(y), dt$L) - x %*% beta
//
//   eta1 = d + (1 - d) * (dt$L > 1e-8)   # eta_{1jl}
//   eta2 = d + (1 - d) * (dt$R < 1e8)    # eta_{2ik}
//   phi = 1/m                             # varphi_j
//
//   gam = matrix(0, p, p)
//   cluster_ids = unique(id)
//   for (i in 1:n) {
//     i_idx = which(id == cluster_ids[i])
//     for (k in i_idx) {
//       Xik = x[k, , drop = FALSE]
//       vik = v[k]
//       ujl = u
//       term1 = matrix(0, 1, p)
//       term1diff = matrix(rep(Xik, N), nrow=N, byrow=FALSE) - x
//       ind1 = (vik < u)
//       w1 = rep(phi, times=m) * eta2[k] / n
//       term1 = colSums(sweep(term1diff , 1, (ind1*w1), FUN='*'))
//
// # second term (denominator + numerator)
//       term2 = matrix(0, 1, p)
//         for (j in 1:n) {
//           j_idx = which(id == cluster_ids[j])
//           for (l in j_idx) {
//             ujl = u[l]
//             ind = as.numeric(vik >= ujl)
//             if(ind == 1) {
//               idx_r = which(v > ujl)
//               if (length(idx_r) > 0) {
//                 mr = 1/m[id[idx_r]]
//                 num = colSums(sweep(matrix(rep(Xik, length(idx_r)), ncol = p, byrow = TRUE) - x[idx_r, , drop = FALSE],
//                                     1, mr, FUN = "*"))
//                 denom = sum(mr)
//                 if (denom > 0) {
//                   term2 = term2 + (phi[j] * eta1[l] / n) * (num / denom)
//                 }
//               }
//             }
//           }
//         }
//         xi_ik = term1 - term2
//         gam = gam + t(xi_ik) %*% xi_ik
//     }
//   }
//   return(gam)
// }


// [[Rcpp::export]]
arma::mat Gamma_func2p_pic_cpp(
    const arma::vec& id,
    const arma::mat& x,
    const arma::vec& L,
    const arma::vec& R,
    const arma::vec& d,
    const arma::vec& beta)
{
  int N = x.n_rows;
  int p = x.n_cols;

  // unique cluster IDs
  arma::vec cluster_ids = arma::unique(id);
  int n = cluster_ids.n_elem;

  // m_j = cluster size
  arma::vec m(n);
  for (int j = 0; j < n; j++) {
    m(j) = sum(id == cluster_ids(j));
  }

  // φ_j = 1/m_j
  arma::vec phi = 1.0 / m;

  // y = max(L, 1e-8)
  arma::vec y = arma::max(L, arma::vec(N).fill(1e-8));
  arma::vec xb = x * beta;

  arma::vec v(N), u_vec(N);

  // define v, u
  for (int i = 0; i < N; i++) {
    if (d(i) == 1) { // exact
      v(i) = std::log(y(i)) - xb(i);
      u_vec(i) = std::log(y(i)) - xb(i);
    } else {
      v(i) = R(i) - xb(i);
      u_vec(i) = L(i) - xb(i);
    }
  }

  // eta1, eta2
  arma::vec eta1(N), eta2(N);
  for (int i = 0; i < N; i++) {
    eta1(i) = d(i) + (1 - d(i)) * (L(i) > 1e-8);
    eta2(i) = d(i) + (1 - d(i)) * (R(i) < 1e8);
  }

  arma::mat gam(p, p, fill::zeros);

  // ======================
  //   Main double loop
  // ======================
  for (int ci = 0; ci < n; ci++) {

    double current_id = cluster_ids(ci);
    arma::uvec i_idx = arma::find(id == current_id);

    // Iterate subject k in cluster i
    for (uword kk = 0; kk < i_idx.n_elem; kk++) {

      int k = i_idx(kk);
      arma::rowvec Xik = x.row(k);
      double vik = v(k);

      arma::uvec ind1_u = (vik < u_vec);

      // w1(j) = (eta2(k)/n) * phi(clust(j))
      arma::vec w1(N);
      double wtemp = eta2(k) / double(n);

      for (int j = 0; j < N; j++) {
        arma::uvec pos = arma::find(cluster_ids == id(j));
        int idx = pos(0); // safe: first match only
        w1(j) = wtemp * phi(idx);
      }

      // term1 = sum_j (ind1(j) * w1(j) * (Xik - x_j))
      arma::mat diff1 = arma::repmat(Xik, N, 1) - x;
      arma::vec  ind1   = arma::conv_to<arma::vec>::from(ind1_u);  // 0/1을 double로
      arma::rowvec term1 = (ind1 % w1).t() * diff1;


      // --------------- term2 ---------------
      arma::rowvec term2(p, fill::zeros);

      for (int cj = 0; cj < n; cj++) {

        arma::uvec j_idx = arma::find(id == cluster_ids(cj));
        double phi_j = phi(cj);

        for (uword ll = 0; ll < j_idx.n_elem; ll++) {

          int l = j_idx(ll);
          double ujl = u_vec(l);

          if (vik >= ujl) {

            arma::uvec idx_r = arma::find(v > ujl);

            if (idx_r.n_elem > 0) {

              arma::vec mr(idx_r.n_elem);
              for (uword rr = 0; rr < idx_r.n_elem; rr++) {

                int r_idx = idx_r(rr);
                arma::uvec pos = arma::find(cluster_ids == id(r_idx));
                int cpos = pos(0);

                mr(rr) = 1.0 / m(cpos);
              }

              arma::rowvec num =
                mr.t() * (arma::repmat(Xik, idx_r.n_elem, 1) - x.rows(idx_r));

              double denom = arma::sum(mr);

              if (denom > 0) {
                term2 += (phi_j * eta1(l) / double(n)) * (num / denom);
              }
            }

          }

        } // l-loop
      } // j cluster

      // xi_ik
      arma::rowvec xi_ik = term1 - term2;

      // accumulate Γ
      gam += xi_ik.t() * xi_ik;
    }
  }

  return gam;
}

// est_funcp = function(ei,ej,xij,si,mi,mj,z=NULL) {
//   n = sqrt(length(si))
//   ind_i = rep(1:n,each=n)
//   if (is.null(z)) z = rep(1,n)
//     zi = z[ind_i]
//   Phi = zi*drop(pnorm((ej-ei)/si))# standard normal distribution function
//   Phi = Phi/(mi*mj)
//     drop(t(xij) %*% Phi)/n
// } # function S(B)

// [[Rcpp::export]]
arma::vec est_funcp_pic_cpp(const arma::vec& ei,
                        const arma::vec& ej,
                        const arma::mat& xij,
                        const arma::vec& si,
                        const arma::vec& mi,
                        const arma::vec& mj,
                        Rcpp::Nullable<arma::vec> z_in = R_NilValue)
{
  int N2 = ei.n_elem;
  int n  = std::sqrt((double)N2);
  int p  = xij.n_cols;

  // -----------------------------
  // 1. z 처리 (Optional)
  // -----------------------------
  arma::vec z;

  if (z_in.isNull()) {
    // z = rep(1, n)
    z = arma::vec(n, arma::fill::ones);
  } else {
    z = Rcpp::as<arma::vec>(z_in);
    if (z.n_elem != n) {
      Rcpp::stop("Length of z must be equal to n.");
    }
  }
  arma::uvec ind_i(N2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ind_i[i*n + j] = i;
    }
  }

  arma::vec zi(N2);
  for (int k = 0; k < N2; ++k) {
    zi[k] = z[ind_i[k]];
  }
  arma::vec Phi(N2);
  for (int k = 0; k < N2; ++k) {
    double ratio = (ej[k] - ei[k]) / si[k];
    Phi[k] = zi[k] * R::pnorm(ratio, 0.0, 1.0, 1, 0);
  }

  Phi = Phi / (mi % mj);

  arma::vec S = (xij.t() * Phi) / (double)n;

  return S;
}


// un_est_funcp = function(ei,ej,xij,si,mi,mj,z=NULL) {
//     n = sqrt(length(si))
//     ind_i = rep(1:n,each=n)
//     if (is.null(z)) z = rep(1,n)
//       zi = z[ind_i]
//     Phi = zi*(ei<ej)  # indicated function
//     Phi = Phi/(mi*mj)
//       drop((t(xij) %*% Phi))/n
//   }

// [[Rcpp::export]]
arma::vec un_est_funcp_pic_cpp(const arma::vec& ei,
                           const arma::vec& ej,
                           const arma::mat& xij,
                           const arma::vec& si,
                           const arma::vec& mi,
                           const arma::vec& mj,
                           const arma::vec& z) {
  int n = std::sqrt(ei.n_elem);

  arma::vec ind = arma::conv_to<arma::vec>::from(ei < ej);

  arma::vec zi = z;
  arma::vec zi_full(ind.n_elem);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      zi_full[i*n + j] = zi[i];
    }
  }
  arma::vec Phi = zi_full % ind;

  arma::vec Phi2 = Phi / (mi % mj);

  arma::mat tmp = xij.each_col() % Phi2;
  arma::vec out = tmp.t() * arma::ones(Phi2.n_elem) / n;

  return out;
}

// Gamma_funcp = function(ei,ej,xij,si,mi,mj,B=200) {
//   n = sqrt(length(si))
//   Shat = t(replicate(B,{
//     un_est_funcp(ei,ej,xij,si,z=rexp(n),mi,mj)*sqrt(n)
//   }))
//
//   cov(Shat)
// }

// [[Rcpp::export]]
arma::mat Gamma_funcp_pic_cpp(const arma::vec& ei,
                          const arma::vec& ej,
                          const arma::mat& xij,
                          const arma::vec& si,
                          const arma::vec& mi,
                          const arma::vec& mj,
                          const int B = 200) {
  int n = std::sqrt(si.n_elem);
  int p = xij.n_cols;

  arma::mat Shat(B, p);

  for (int b = 0; b < B; b++) {

    arma::vec z(n);
    for (int i = 0; i < n; i++) {
      z[i] = R::rexp(1.0);
    }

    arma::vec score = un_est_funcp_pic_cpp(ei, ej, xij, si, mi, mj, z);
    Shat.row(b) = (score.t() * std::sqrt((double)n));
  }
  arma::rowvec mean_row = arma::mean(Shat,0);
  arma::mat centered = Shat.each_row() - mean_row;
  arma::mat Gamma = (centered.t() * centered ) / (B-1) ;

  return Gamma;
}

// up_Amatp = function(ei,ej,xij,si,mi,mj,n) {
//   phi = drop(dnorm((ej-ei)/si))
//   phi = phi/(mi*mj)
//   crossprod(xij*phi/si,xij)/n
// } # hessian

// [[Rcpp::export]]
arma::mat up_Amatp_pic_cpp(const arma::vec& ei,
                      const arma::vec& ej,
                      const arma::mat& xij,
                      const arma::vec& si,
                      const arma::vec& mi,
                      const arma::vec& mj) {
  int n = std::sqrt(ei.n_elem);
  arma::vec ratio = (ej - ei) / si;
  arma::vec phi(ratio.n_elem);

  for (uword k = 0; k < ratio.n_elem; ++k) {
    phi[k] = R::dnorm4(ratio[k], 0.0, 1.0, 0);
  }
  arma::vec mimj = mi % mj;
  arma::vec phi2 = phi /mimj;
  arma::vec weight = phi2 / si;
  arma::mat A = (xij.each_col() % weight).t() *  xij/n;

  return A;
}



// Sigma_funcp = function(ei,ej,xij,Gamma,si,mi,mj,n) {
//   A = up_Amatp(ei,ej,xij,si,mi,mj,n)
//   min_eig = min(eigen(A)$value)
//   if(min_eig < 1e-10) {
//     Ainv = ginv(A + diag(1e-05, p,p))
//   } else {
//     Ainv = solve(A)
//   }
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_funcp_pic_cpp(const arma::vec& ei,
                         const arma::vec& ej,
                         const arma::mat& xij,
                         const arma::mat& Gamma,
                         const arma::vec& si,
                         const arma::vec& mi,
                         const arma::vec& mj) {
  int n = std::sqrt(ei.n_elem);
  arma::mat A = up_Amatp_pic_cpp(ei, ej, xij, si, mi, mj);
  arma::mat Ainv;
  bool ok = arma::inv_sympd(Ainv, A);
  if(!ok) {
    Ainv = arma::inv(A);
  }
  arma::mat Sigma = Ainv * Gamma * Ainv / n;
  return Sigma;
}

// Sigma_funcp_pen = function(ei,ej,xij,Gamma,si,mi,mj,n, tilde_beta, lambda) {
//   A = up_Amatp(ei,ej,xij,si,mi,mj,n) + diag((2*lambda) / (tilde_beta^2 + 1e-05))
//   min_eig = min(eigen(A)$value)
//   if(min_eig < 1e-10) {
//     Ainv = ginv(A + diag(1e-05, p,p))
//   } else {
//     Ainv = solve(A)
//   }
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_funcp_pen_pic_cpp(const arma::vec&ei,
                              const arma::vec&ej,
                              const arma::mat& xij,
                              const arma::mat& Gamma,
                              const arma::vec& si,
                              const arma::vec& mi,
                              const arma::vec& mj,
                              const int n,
                              const arma::vec& tilde_beta,
                              const double lambda) {
  int p = tilde_beta.n_elem;
  arma::mat A = up_Amatp_pic_cpp(ei, ej, xij, si, mi, mj);
  arma::vec pen = (2.0 * lambda) / (arma::square(tilde_beta) + 1e-5) ;
  A += arma::diagmat(pen);

  arma::vec eigvals; // confirm eigen value
  arma::eig_sym(eigvals, A);
  double min_eig = eigvals.min();

  arma::mat Ainv;

  if (min_eig < 1e-10) {
    arma::mat A_eps = A + 1e-5 * arma::eye(p,p);
    Ainv = arma::pinv(A_eps);
  } else {
    Ainv = arma::inv(A);
  }

  arma::mat out = (Ainv * Gamma * Ainv) / n;
  return out;
}

// sifunc = function(xij, Sigma) {
//   si = pmax(0.0001, rowSums((xij %*% Sigma) * xij))
//
//   sqrt(si)
// } # sifunc = ((Xi-Xj)*sigma*(Xi-Xj))^(1/ 2)

// [[Rcpp::export]]
arma::vec sifunc_pic_cpp(const arma::mat& xij,
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

// is_aftp = function(dt, se, beta = NULL, maxit = 10, tol = 1e-2) {
//   y = pmax(dt$L, 1e-08)
//   d = dt$delta
//   x = as.matrix(dt[, -c(1:4)])
//   n = nrow(dt)
//   id = dt$id
//   p = ncol(x)
//   m = unname(table(id))
//
//   L = log(y); R = log(dt$R)
//     Sigma = diag(p) / n
//
//   if (is.null(beta)) {
//     beta = lm(log(y) ~ x, weights = d)$coef[-1]
//     iterative = TRUE
//   } else {
//     iterative = FALSE
//   }
//
//   id_i = rep(1:n, each = n)
//     id_j = rep(1:n, times = n)
//     m = rep(m, m); mi = m[id_i]; mj = m[id_j]
//     xi = x[id_i, ]; xj = x[id_j, ]
//     xij = xi - xj
//
//     yi = ifelse(d == 1, log(y), R)[id_i]
//     yj = ifelse(d == 1, log(y), L)[id_j]
//
//     ei = yi - xi %*% beta
//     ej = yj - xj %*% beta
//     si = sifunc(xij, Sigma)
//
//       iter = 0; err = 10
//
//       while (iter < maxit & err > tol & iterative) {
//         A = up_Amatp(ei, ej, xij, si, mi, mj, n)
//         S = est_funcp(ei, ej, xij, si, mi, mj)
//         new_beta = drop(beta - S %*% solve(A))
//         err = norm(new_beta - beta, "2")
//         beta = new_beta
//         ei = yi - xi %*% beta
//         ej = yj - xj %*% beta
//         if (se == "ZL") Gamma = Gamma_funcp(ei, ej, xij, si, mi, mj)
//           if (se == "CF") Gamma = Gamma_func2p(dt, beta)
//             Sigma = Sigma_funcp(ei, ej, xij, Gamma, si, mi, mj, n)
//             si = sifunc(xij, Sigma)
//
//             iter = iter + 1
//       }
//
//       A = up_Amatp(ei, ej, xij, si, mi, mj, n)
//         S = est_funcp(ei, ej, xij, si, mi, mj)
//
//         obj = list(beta = new_beta,
//                    hess = A,
//                    grad = S)
//         return(obj)
//
// }

// [[Rcpp::export]]
Rcpp::List is_aftp_pic_cpp(Rcpp::DataFrame dt,
                       std::string se,
                       Rcpp::Nullable<arma::vec> beta_in = R_NilValue,
                       int maxit = 10,
                       double tol = 1e-2)
{
  using arma::vec;
  using arma::mat;

  // 기본 세팅 ----------------------------------------------------
  int n    = dt.nrows();
  int ncol = dt.size();
  int p    = ncol - 4;  // L, R, delta, id 제외

  vec L    = Rcpp::as<vec>(dt["L"]);
  vec Rvec = Rcpp::as<vec>(dt["R"]);
  vec d    = Rcpp::as<vec>(dt["delta"]);
  vec id   = Rcpp::as<vec>(dt["id"]);

  // design matrix x: 5번째 컬럼부터 covariate
  mat x(n, p);
  for (int j = 0; j < p; ++j) {
    x.col(j) = Rcpp::as<vec>(dt[j + 4]);
  }

  vec y    = arma::max(L, vec(n).fill(1e-8));
  vec logy = arma::log(y);

  mat Sigma = arma::eye<mat>(p, p) / static_cast<double>(n);

  // 초기값 beta ---------------------------------------------------
  vec beta;
  bool iterative;
  if (beta_in.isNull()) {
    // weighted LS 초기값 (right-censoring에서 d를 weight로)
    vec w     = d;
    vec sqrtw = arma::sqrt(w);

    mat Xw = x.each_col() % sqrtw;
    vec yw = logy % sqrtw;

    beta      = arma::solve(Xw.t() * Xw, Xw.t() * yw);
    iterative = true;
  } else {
    beta      = Rcpp::as<vec>(beta_in);
    iterative = false;
  }

  // cluster size m_per_sub ---------------------------------------
  vec m_per_sub(n);
  {
    int i = 0;
    while (i < n) {
      double cur_id = id[i];
      int start = i;
      while (i < n && id[i] == cur_id) {
        i++;
      }
      int m_cluster = i - start;
      for (int k = start; k < i; ++k) {
        m_per_sub[k] = m_cluster;
      }
    }
  }

  // pair index ---------------------------------------------------
  int N2 = n * n;
  arma::uvec id_i(N2), id_j(N2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int idx = i * n + j;
      id_i[idx] = i;
      id_j[idx] = j;
    }
  }
  vec mi = m_per_sub.elem(id_i);
  vec mj = m_per_sub.elem(id_j);

  // x_{ij} = x_i - x_j -------------------------------------------
  mat xi(N2, p), xj(N2, p);
  for (int k = 0; k < N2; ++k) {
    xi.row(k) = x.row(id_i[k]);
    xj.row(k) = x.row(id_j[k]);
  }
  mat xij = xi - xj;

  // y_i, y_j (partly interval-censored rule) ---------------------
  vec yi(N2), yj(N2);
  for (int k = 0; k < N2; ++k) {
    int ii = id_i[k];
    int jj = id_j[k];

    double yi_val = (d[ii] == 1.0) ? logy[ii] : std::log(Rvec[ii]); // event vs right / interval
    double yj_val = (d[jj] == 1.0) ? logy[jj] : std::log(L[jj]);    // 여기는 네가 쓰던 rule 그대로 둠

    yi[k] = yi_val;
    yj[k] = yj_val;
  }

  vec ei = yi - (xi * beta);
  vec ej = yj - (xj * beta);
  vec si = sifunc_pic_cpp(xij, Sigma);

  // 반복 알고리즘 ------------------------------------------------
  int iter   = 0;
  double err = 10.0;

  mat A(p, p);
  vec S(p);
  mat Gamma(p, p, arma::fill::zeros);

  while (iter < maxit && err > tol && iterative) {

    // 1. Hessian, Score
    A = up_Amatp_pic_cpp(ei, ej, xij, si, mi, mj);
    S = est_funcp_pic_cpp(ei, ej, xij, si, mi, mj);

    // 2. beta 업데이트
    vec new_beta = beta - arma::solve(A, S);
    err  = arma::norm(new_beta - beta, 2);
    beta = new_beta;

    // 3. residual 업데이트
    ei = yi - (xi * beta);
    ej = yj - (xj * beta);

    // 4. Gamma (variance of score) 선택
    if (se == "ZL") {
      Gamma = Gamma_funcp_pic_cpp(ei, ej, xij, si, mi, mj);
    } else if (se == "CF") {
      Gamma = Gamma_func2p_pic_cpp(id, x, L, Rvec, d, beta);
    }

    // 5. Sandwich Sigma
    Sigma = Sigma_funcp_pic_cpp(ei, ej, xij, Gamma, si, mi, mj);
    si    = sifunc_pic_cpp(xij, Sigma);

    iter++;
  }

  // 마지막 Hessian, Score 한 번 더 계산 -------------------------
  A = up_Amatp_pic_cpp(ei, ej, xij, si, mi, mj);
  S = est_funcp_pic_cpp(ei, ej, xij, si, mi, mj);

  return Rcpp::List::create(
    Rcpp::Named("beta")  = beta,
    Rcpp::Named("hess")  = A,
    Rcpp::Named("grad")  = S,
    Rcpp::Named("Sigma") = Sigma,
    Rcpp::Named("iter")  = iter
  );
}



// est_func = function(ei,ej,xij,eta2i, eta1j ,si,z=NULL) {
//   n = sqrt(length(si))
//   ind_i = rep(1:n,each=n)
//   if (is.null(z)) z = rep(1,n)
//     zi = z[ind_i]
//   Phi = zi*drop(pnorm((ej-ei)/si)) # standard normal distribution function
//   d = eta2i*eta1j
//   drop(t(d*xij) %*% Phi)/n
// } # first derivative func

// [[Rcpp::export]]
arma::vec est_func_pic_cpp(const arma::vec& ei,
                       const arma::vec& ej,
                       const arma::mat& xij,
                       const arma::vec& eta2i,
                       const arma::vec& eta1j,
                       const arma::vec& si) {
  int N2 = ei.n_elem;
  int p = xij.n_cols;
  int n = std::sqrt((double)N2);

  arma::vec Phi(N2);
  for (int k = 0; k < N2; k++)
    Phi[k] = R::pnorm((ej[k] - ei[k]) / si[k], 0, 1, 1, 0);

  arma::vec d = eta2i % eta1j;

  return (xij.each_col() % (d % Phi / n)).t() * arma::ones(N2);
}

// un_est_func = function(ei,ej,xij,eta2i,eta1j,si,z=NULL) {
//   n = sqrt(length(si))
//   ind_i = rep(1:n,each=n)
//   if (is.null(z)) z = rep(1,n)
//     zi = z[ind_i]
//   Phi = zi*(ei<ej) # indicated function
//   d = eta2i*eta1j
//   drop(t(d*xij) %*% Phi)/n
// }

// [[Rcpp::export]]
arma::vec un_est_func_pic_cpp(const arma::vec& ei,
                           const arma::vec& ej,
                           const arma::mat& xij,   // (n^2) x p
                           const arma::vec& eta2i,
                           const arma::vec& eta1j,
                           const arma::vec& si,    // length n^2 (only for n)
                           const arma::vec& z) {   // length n

  const arma::uword N = ei.n_elem;
  const arma::uword n = (arma::uword) std::llround(std::sqrt((double)N));

  // zi_full[k] = z[i] where i = k/n  (rep(1:n, each=n))
  vec zi_full(N);
  for (arma::uword k = 0; k < N; ++k) zi_full[k] = z[k / n];

  // Phi = zi_full * (ei < ej)
  vec Phi = zi_full % arma::conv_to<vec>::from(ei < ej);

  // d = eta2i * eta1j
  vec w = (eta2i % eta1j) % Phi;

  // drop(t(d*xij) %*% Phi)/n  ==  (xij.t() * w)/n
  return (xij.t() * w) / (double)n;
}

// up_Amat = function(ei,ej,xij,eta2i,eta1j,si) {
//   n = sqrt(length(si))
//   phi = drop(dnorm((ej-ei)/si))
//   d = eta2i * eta1j
//   crossprod(d*xij*phi/si,xij)/n
// } # hessian

// Gamma_func = function(ei,ej,xij,eta2i,eta1j,si,B=200) {
//   n = sqrt(length(si))
//   Shat = t(replicate(B,{
//     un_est_func(ei,ej,xij,eta2i,eta1j,si,z=rexp(n))*sqrt(n)
//   }))
//
//   cov(Shat)
// }

// [[Rcpp::export]]
arma::mat up_Amat_pic_cpp(const arma::vec& ei,
                      const arma::vec& ej,
                      const arma::mat& xij,
                      const arma::vec& eta2i,
                      const arma::vec& eta1j,
                      const arma::vec& si) {
  const int N2 = ei.n_elem;
  const int p  = xij.n_cols;
  const int n  = static_cast<int>(std::sqrt(static_cast<double>(N2)));

  arma::vec phi(N2);

  for (int k = 0; k < N2; k++) {
    const double s = si(k);
    const double z = (ej(k) - ei(k)) / s;
    phi(k) = R::dnorm(z, 0.0, 1.0, 0);
  }

  arma::vec d = eta2i % eta1j;
  arma::vec w = (d % phi) / si;
  arma::mat W = xij.each_col() % w;

  return (W.t() * xij) / static_cast<double>(n);
}

// [[Rcpp::export]]
arma::mat Gamma_func_pic_cpp(const arma::vec& ei,
                         const arma::vec& ej,
                         const arma::mat& xij,
                         const arma::vec& eta2i,
                         const arma::vec& eta1j,
                         const arma::vec& si,
                         int B = 200) {
  int N2 = ei.n_elem;
  int p = xij.n_cols;
  int n = std::sqrt((double)N2);

  arma::mat Shat(B, p);

  for (int b = 0; b < B; b++) {
    arma::vec z(n);
    for (int i = 0; i < n; i++) z[i] = R::rexp(1.0);

    arma::vec score = un_est_func_pic_cpp(ei, ej, xij, eta2i, eta1j, si, z);
    Shat.row(b) = (score.t() * std::sqrt((double)n));
  }
  arma::rowvec mean_row = arma::mean(Shat, 0);
  arma::mat centered = Shat.each_row() - mean_row;

  return (centered.t() * centered) / (B-1);
}


// Gamma_func2 = function(dt, beta) {
//   n = length(unique(dt$id))
//   id = dt$id
//   m = as.numeric(table(id))
//   N = sum(m)
//   x = as.matrix(dt[, -(1:4)])
//   p = ncol(x)
//   y = pmax(dt$L, 1e-8)
//   d = dt$delta
//   v = ifelse(d == 1, log(y), dt$R) - x %*% beta
//   u = ifelse(d == 1, log(y), dt$L) - x %*% beta
//
//   eta1 = d + (1 - d) * (dt$L > 1e-8)   # eta_{1jl}
//   eta2 = d + (1 - d) * (dt$R < 1e8)    # eta_{2ik}
//
//   gam = matrix(0, p, p)
//   cluster_ids = unique(id)
//   for (i in 1:n) {
//     i_idx = which(id == cluster_ids[i])
//     for (k in i_idx) {
//       Xik = x[k, , drop = FALSE]
//       vik = v[k]
//       ujl = u
//       term1 = matrix(0, 1, p)
//       term1diff = matrix(rep(Xik, N), nrow=N, byrow=FALSE) - x
//       ind1 = (vik < u)
//       w1 = eta2[k] / n
//       term1 = colSums(sweep(term1diff , 1, (ind1*w1), FUN='*'))
//
// # second term (denominator + numerator)
//       term2 = matrix(0, 1, p)
//         for (j in 1:n) {
//           j_idx = which(id == cluster_ids[j])
//           for (l in j_idx) {
//             ujl = u[l]
//             ind = as.numeric(vik >= ujl)
//             if(ind == 1) {
//               idx_r = which(v > ujl)
//               if (length(idx_r) > 0) {
//                 num = colSums(matrix(rep(Xik, length(idx_r)), ncol = p, byrow = TRUE) -
//                   x[idx_r, , drop = FALSE])
//                 denom = length(idx_r)
//                 if (denom > 0) {
//                   term2 = term2 + (eta1[l] / n) * (num / denom)
//                 }
//               }
//             }
//           }
//         }
//         xi_ik = term1 - term2
//         gam = gam + t(xi_ik) %*% xi_ik
//     }
//   }
//   return(gam)
// }

// [[Rcpp::export]]
arma::mat Gamma_func2_pic_cpp(Rcpp::DataFrame dt,
                           const arma::vec& beta) {
  arma::vec id = dt["id"];
  arma::vec L = dt["L"];
  arma::vec Rvec = dt["R"];
  arma::vec d = dt["delta"];

  int N = id.n_elem;

  int p = dt.size() - 4;
  arma::mat x(N, p);
  for (int j = 0; j < p; j++) {
    x.col(j) = Rcpp::as<arma::vec>(dt[j + 4]);
  }

  arma::vec uniq_id = arma::unique(id);
  int n = uniq_id.n_elem;

  arma::vec m_per_sub(N);
  {
    int i = 0;
    while (i < N) {
      const double cur = id(i);
      const int start = i;
      while (i < N && id(i) == cur) i++;
      const int m = i - start;
      for (int k = start; k < i; k++) m_per_sub[k] = m;
    }
  }

  arma::vec eta1(N), eta2(N);
  arma::vec y = arma::max(L, arma::vec(N).fill(1e-8));
  arma::vec v(N), u(N);
  for (int i = 0; i < N; i++) {
    eta1(i) = d(i) + (1.0 - d(i)) * (L(i) > 1e-8 ? 1.0 : 0.0);
    eta2(i) = d(i) + (1.0 - d(i)) * (Rvec(i) < 1e8 ? 1.0 : 0.0);
    double logy_i = std::log(y[i]);
    if (d[i] == 1) {
      v[i] = logy_i - arma::dot(x.row(i), beta);
      u[i] = logy_i - arma::dot(x.row(i), beta);
    } else {
      v[i] = Rvec[i] - arma::dot(x.row(i), beta);
      u[i] = L[i] - arma::dot(x.row(i), beta);
      }
  }
  arma::mat gam(p, p, arma::fill::zeros);

  for (int ii = 0; ii < n ; ii++) {
    double cid = uniq_id(ii);
    std::vector<int> i_idx;
    for (int k = 0; k < N; k++) if (id[k] == cid) i_idx.push_back(k);

    for (int k : i_idx) {
      arma::rowvec Xik = x.row(k);
      double vik = v[k];

      arma::mat diff = arma::repmat(Xik, N, 1) - x;

      arma::vec ind1(N);
      for (int j = 0; j < N; j++)
        ind1(j) = (vik < u(j)) ? 1.0 : 0.0;

      double w1 = eta2[k] / n;
      arma::rowvec term1 = (ind1.t() * diff) * w1;

      arma::rowvec term2(p, arma::fill::zeros);

      for (int jj = 0; jj < n; jj++) {
        double cid2= uniq_id[jj];
        std::vector<int> j_idx;
        for (int t = 0; t < N; t++)
          if (id[t] == cid2) j_idx.push_back(t);

        for (int l : j_idx) {
          double ujl = u[l];

          if (vik >= ujl) {
            std::vector<int> idx_r;
            for (int r = 0; r < N; r++)
              if (v[r] > ujl) idx_r.push_back(r);

            int denom = idx_r.size();
            if (denom > 0) {
              arma::rowvec num(p, arma::fill::zeros);
              for (int r : idx_r)
                num += (Xik - x.row(r));

              num /= denom;

              term2 += (eta1[l] / n) * num;
            }
          }
        }
      }
      arma::rowvec xi_ik = term1 - term2;
      gam += xi_ik.t() * xi_ik;
    }
  }
  return gam;
}


// Sigma_func = function(ei,ej,xij,eta2i,eta1j,Gamma,si) {
//   n = sqrt(length(si))
//   A = up_Amat(ei,ej,xij,eta2i,eta1j,si)
//   min_eig = min(eigen(A)$value)
//   if(min_eig < 1e-10) {
//     Ainv = ginv(A + diag(1e-05, p,p))
//   } else {
//     Ainv = solve(A)
//   }
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_func_pic_cpp(const arma::vec& ei,
                         const arma::vec& ej,
                         const arma::mat& xij,
                         const arma::vec& eta2i,
                         const arma::vec& eta1j,
                         const arma::mat& Gamma,
                         const arma::vec& si) {
  int p = Gamma.n_rows;
  int N2 = ei.n_elem;
  int n = std::sqrt((double)N2);

  arma::mat A = up_Amat_pic_cpp(ei, ej, xij, eta2i, eta1j, si);

  arma::vec eval;
  arma::eig_sym(eval, A);

  arma::mat Ainv;
  if (eval.min() < 1e-10) {
    Ainv = arma::pinv(A + 1e-5 * arma::eye(p,p));
  } else {
    Ainv = arma::inv(A);
  }
  return Ainv * Gamma * Ainv / n;
}


// Sigma_func_pen = function(ei,ej,xij,eta2i,eta1j,Gamma,si, tilde_beta, lambda) {
//   n = sqrt(length(si))
//   A = up_Amat(ei,ej,xij,eta2i,eta1j,si) + diag((2*lambda) / (tilde_beta^2 + 1e-05))
//   min_eig = min(eigen(A)$value)
//   if(min_eig < 1e-10) {
//     Ainv = ginv(A + diag(1e-05, p,p))
//   } else {
//     Ainv = solve(A)
//   }
//   Ainv %*% Gamma %*% Ainv/n
// }

// [[Rcpp::export]]
arma::mat Sigma_func_pen_pic_cpp(const arma::vec& ei,
                             const arma::vec& ej,
                             const arma::mat& xij,
                             const arma::vec& eta2i,
                             const arma::vec& eta1j,
                             const arma::mat& Gamma,
                             const arma::vec& si,
                             const arma::vec& tilde_beta,
                             const double lambda) {
  int p = tilde_beta.n_elem;
  int N2 = si.n_elem;
  int n = std::sqrt((double)N2);

  arma::mat A = up_Amat_pic_cpp(ei, ej, xij, eta2i, eta1j, si);
  arma::vec pen = (2.0 * lambda) / (arma::square(tilde_beta) + 1e-5);
  A += arma::diagmat(pen);

  arma::vec eigvals;
  arma::eig_sym(eigvals, A);
  double min_eig = eigvals.min();

  arma::mat Ainv;

  if (min_eig < 1e-10) {
    arma::mat A_eps = A + 1e-5 * arma::eye(p,p);
    Ainv = arma::pinv(A_eps);
  } else {
    Ainv = arma::inv(A);
  }

  arma::mat out = Ainv * Gamma * Ainv / n;
  return out;
}

//
// is_aft = function(dt,se) {
//   y = pmax(dt$L, 1e-08)
//   d = dt$delta
//   x = as.matrix(dt[, -c(1:4)])
//   n = nrow(dt)
//   p = ncol(x)
//   L = log(y); R = log(dt$R)
//     old_beta = lm(log(y) ~ x, weights = d)$coef[-1]
//   Sigma = diag(p) / n
//   id_i = rep(1:n, each = n)
//     id_j = rep(1:n, times = n)
//     xi = x[id_i, ]; xj = x[id_j, ]
//     xij = xi - xj
//     v = ifelse(d == 1, log(y), R) - x %*% old_beta ; u = ifelse(d == 1, log(y), L) - x %*% old_beta
//     yi = ifelse(d == 1, log(y), R)[id_i] ;yj = ifelse(d == 1, log(y), L)[id_j]
//     ei = yi - xi %*% old_beta ; ej = yj - xj %*% old_beta
//     eta1 = d + (1 - d) * (dt$L > 1e-8) ; eta2 = d + (1 - d) * (dt$R < 1e8)
//       eta2i = eta2[id_i] ; eta1j = eta1[id_j]
//       si = sifunc(xij, Sigma)
//
//         maxit = 100; iter = 0; err = 10; tol = 1e-2
//         while (iter < maxit & tol < err) {
//           A = up_Amat(ei,ej,xij,eta2i,eta1j,si)
//           S = est_func(ei,ej,xij,eta2i,eta1j,si)
//           new_beta = drop(old_beta - S %*% solve(A)) # newton raphson method
//           err = norm(new_beta - old_beta,"2")
// # check for convergence of beta
//           old_beta = new_beta
//           ei = yi - xi%*%old_beta
//           ej = yj - xj%*%old_beta
//           v = ifelse(d == 1, log(y), R) - x %*% old_beta
//           u = ifelse(d == 1, log(y), L) - x %*% old_beta
//           if (se == "ZL") Gamma = Gamma_func(ei,ej,xij,eta2i,eta1j,si)
//             if (se == "CF") Gamma = Gamma_func2(dt, old_beta)
//               Sigma = Sigma_func(ei,ej,xij,eta2i,eta1j,Gamma,si)
//               si = sifunc(xij,Sigma)
//               iter = iter + 1
//         }
//         obj = list(beta = new_beta,
//                    hess = A,
//                    grad = S)
//           return(obj)
// }

// [[Rcpp::export]]
Rcpp::List is_aft_pic_cpp(Rcpp::DataFrame dt, std::string se) {
  arma::vec L = dt["L"];
  arma::vec Rvec = dt["R"];
  arma::vec d = dt["delta"];

  int n = L.n_elem;

  int p = dt.size() - 4;
  arma::mat x(n, p);
  for (int j = 0; j < p; j++) {
    x.col(j) = Rcpp::as<arma::vec>(dt[j + 4]);
  }


  arma::vec y = arma::max(L, arma::vec(n).fill(1e-8));
  arma::vec logy = arma::log(y);

  arma::mat Xw = x.each_col() % arma::sqrt(d);
  arma::vec yw = logy % arma::sqrt(d);
  arma::vec beta = arma::solve(Xw.t() * Xw, Xw.t() * yw);

  arma::mat Sigma = arma::eye(p,p) / n;

  int N2 = n * n;
  arma::uvec id_i(N2), id_j(N2);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      int k = i * n + j;
      id_i[k] = i;
      id_j[k] = j;
    }
  arma::mat xi(N2, p), xj(N2, p);
  for (int k = 0; k < N2; k++) {
    xi.row(k) = x.row(id_i[k]);
    xj.row(k) = x.row(id_j[k]);
  }
  arma::mat xij = xi - xj;

  arma::vec yi(N2), yj(N2);
  for (int k = 0; k < N2; k++) {
    int ii = id_i[k];
    int jj = id_j[k];
    yi[k] = (d[ii] == 1 ? logy[ii] : std::log(Rvec[ii]));
    yj[k] = (d[jj] == 1 ? logy[jj] : std::log(L[jj]));
  }
  arma::vec eta1 = d + (1-d) % (L > 1e-8);
  arma::vec eta2 = d + (1-d) % (Rvec < 1e8);
  arma::vec eta1j = eta1.elem(id_j);
  arma::vec eta2i = eta2.elem(id_i);

  arma::vec ei = yi - (xi * beta);
  arma::vec ej = yj - (xj * beta);

  arma::vec si = sifunc_pic_cpp(xij, Sigma);

  int iter=0;
  double err = 10, tol = 1e-2;
  arma::mat A;
  arma::vec S;
  arma::mat Gamma;

  while (iter < 100 && err > tol) {
    A = up_Amat_pic_cpp(ei, ej, xij, eta2i, eta1j, si);
    S = est_func_pic_cpp(ei, ej, xij, eta2i, eta1j, si);
    arma::vec new_beta = beta - arma::solve(A,S);
    err = arma::norm(new_beta - beta, 2);
    beta = new_beta;

    ei = yi - (xi * beta);
    ej = yj - (xj * beta);

    if (se == "ZL")
      Gamma = Gamma_func_pic_cpp(ei, ej, xij, eta2i, eta1j, si);
    else if (se == "CF")
      Gamma = Gamma_func2_pic_cpp(dt, beta);

    Sigma = Sigma_func_pic_cpp(ei, ej, xij, eta2i, eta1j, Gamma, si);
    si = sifunc_pic_cpp(xij, Sigma);
    iter++;
  }
  A = up_Amat_pic_cpp(ei, ej, xij, eta2i, eta1j, si);
  S = est_func_pic_cpp(ei, ej, xij, eta2i, eta1j, si);

  return Rcpp::List::create(
    _["beta"] = beta,
    _["hess"] = A,
    _["grad"] = S
  );

}
