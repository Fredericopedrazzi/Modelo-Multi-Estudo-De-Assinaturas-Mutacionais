// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// ----------------- helpers -----------------
inline bool is_fixed(const std::vector<int>& fix_cols, int one_based_col) {
    return std::find(fix_cols.begin(), fix_cols.end(), one_based_col) != fix_cols.end();
}
inline double lgamma_R(double x) { return R::lgammafn(x); }

// ----------------- update_A (MH/temperado) -----------------
std::vector<arma::mat> update_A(
    const arma::mat& P,
    const std::vector<arma::mat>& E_list,
    std::vector<arma::mat> A_list,
    const std::vector<arma::imat>& M_list,
    const std::vector<arma::imat>& W_list,
    int S,
    int n_signatures,
    double l,
    const std::vector<int>& fix_colsA
) {
    for (int s = 0; s < S; s++) {

        arma::mat As = A_list[s];
        const arma::mat& Es = E_list[s];
        const arma::imat& Ms = M_list[s];
        const arma::imat& Ws = W_list[s];

        // lgamma(M + 1)
        arma::mat Ms_d = arma::conv_to<arma::mat>::from(Ms);
        arma::mat lgamma_term = arma::lgamma(Ms_d + 1.0);

        // para equivalência exata, não usar permutação
        for (int m = 0; m < n_signatures; m++) {

            if (is_fixed(fix_colsA, m + 1)) continue;

            double probl = R::rbeta(l, 3.0 * l);

            arma::mat A1 = As;   // exatamente como no R
            arma::mat A0 = As;
            A0(m, m) = probl;

            // (P %*% A) %*% E   — equivalente ao (P*(A*E))
            arma::mat P_A1 = (P * A1) * Ws(1,1);
            arma::mat P_A0 = (P * A0) * Ws(1,1);

            arma::mat X1 = P_A1 * Es;
            arma::mat X0 = P_A0 * Es;

            // ifelse(is.infinite(log(X)), 0, log(X))
            arma::mat log_X1 = arma::log(X1);
            arma::mat log_X0 = arma::log(X0);

            log_X1.elem(arma::find_nonfinite(log_X1)).zeros();
            log_X0.elem(arma::find_nonfinite(log_X0)).zeros();

            // prob1 = sum( -(X1) + M*log(X1) - lgamma(M+1) )
            double prob1 =
                arma::accu( -X1 + Ms_d % log_X1 - lgamma_term );

            double prob2 =
                arma::accu( -X0 + Ms_d % log_X0 - lgamma_term );

            // razao = exp(prob2 - prob1)
            double razao = std::exp(prob2 - prob1);

            if (R::runif(0.0, 1.0) < razao) {
                As(m, m) = probl;
            }
        }

        A_list[s] = As;
    }

    return A_list;
}

// ----------------- log-sum-exp helper -----------------
inline double log_sum_exp(const std::vector<double> &log_probs) {
    double max_log = *std::max_element(log_probs.begin(), log_probs.end());
    double sum_exp = 0.0;
    for (double v : log_probs) sum_exp += std::exp(v - max_log);
    return max_log + std::log(sum_exp);
}

// ----------------- update_A (versão com log-sum-exp) -----------------
std::vector<arma::mat> update_A_2(const arma::mat& P,
                                const std::vector<arma::mat>& E_list,
                                std::vector<arma::mat> A_list,
                                const std::vector<arma::imat>& M_list,
                                const std::vector<arma::imat>& W_list,
                                int S,
                                int n_signatures,
                                double l,
                                const std::vector<int>& fix_colsA) {

    for (int s = 0; s < S; s++) {
        arma::mat As = A_list[s];
        const arma::mat& Es = E_list[s];
        const arma::imat& Ms = M_list[s];
        const arma::imat& Ws = W_list[s];

        // soma de lgamma, como no R
        double lgamma_sum = 0.0;
        for (uword r = 0; r < Ms.n_rows; r++)
            for (uword c = 0; c < Ms.n_cols; c++)
                lgamma_sum += lgamma_R(static_cast<double>(Ms(r, c)) + 1.0);

        for (int m = 0; m < n_signatures; m++) {
            if (is_fixed(fix_colsA, m + 1)) continue;

            arma::mat A1 = As;
            arma::mat A0 = As;
            A1(m, m) = 1.0;
            A0(m, m) = 0.0;

            arma::mat X1 = P * (A1 * Es)* Ws(1, 1);   //* Ws(1, 1)
            arma::mat X0 = P * (A0 * Es)* Ws(1, 1);

            arma::mat Ms_d = arma::conv_to<arma::mat>::from(Ms);

            // log-probabilidades (iguais ao R)
            double log_prob1 = std::log(0.25) + arma::accu(-X1 + Ms_d % arma::log(X1)) - lgamma_sum;
            double log_prob2 = std::log(0.75) + arma::accu(-X0 + Ms_d % arma::log(X0)) - lgamma_sum;

            
            // log-sum-exp para normalizar
            std::vector<double> log_probs = {log_prob1, log_prob2};
            double denom = log_sum_exp(log_probs);
            double prob = std::exp(log_prob1 - denom);

            // limita probabilidade
            prob = std::min(1.0, std::max(0.0, prob));

            // amostra binomial
            As(m, m) = R::rbinom(1, prob);
        }

        A_list[s] = As;
    }

    return A_list;
}




// ----------------- update_Fi_Z -----------------
void update_Fi_Z(const arma::mat& P,
                 const std::vector<arma::mat>& A_list,
                 const std::vector<arma::mat>& E_list,
                 std::vector<arma::cube>& Z_list,
                 std::vector<arma::cube>& Fi_list,
                 const std::vector<arma::imat>& M_list,
                 int S,
                 int n_signatures,
                 const std::vector<int>& j_list) {

    for (int s = 0; s < S; s++) {
        const arma::mat& As = A_list[s];
        const arma::mat& Es = E_list[s];
        arma::cube& Zs = Z_list[s];
        arma::cube& Fis = Fi_list[s];
        const arma::imat& Ms = M_list[s];

        arma::mat PE = P * (As * Es);

        // Fi
        for (int m = 0; m < n_signatures; m++) {
            double Amm = As(m, m);
            arma::vec Pm = P.col(m);
            arma::rowvec Em = Es.row(m);
            Fis.slice(m) = (Pm * (Amm * Em)) / PE;
        }

        // Z ~ Multinomial
        for (uword k = 0; k < Zs.n_rows; k++) {
            for (uword g = 0; g < Zs.n_cols; g++) {
                int size = Ms(k, g);
                if (size <= 0) {
                    for (int m = 0; m < n_signatures; m++) Zs(k, g, m) = 0.0;
                    continue;
                }
                NumericVector prob(n_signatures);
                double sumprob = 0.0;
                for (int m = 0; m < n_signatures; m++) {
                    double p = Fis(k, g, m);
                    if (!R_finite(p) || p < 0.0) p = 0.0;
                    prob[m] = p;
                    sumprob += p;
                }
                if (sumprob <= 0.0) {
                    for (int m = 0; m < n_signatures; m++) Zs(k, g, m) = 0.0;
                    continue;
                }
                for (int m = 0; m < n_signatures; m++) prob[m] /= sumprob;

                IntegerVector draw(n_signatures);
                R::rmultinom(size, prob.begin(), n_signatures, INTEGER(draw));
                for (int m = 0; m < n_signatures; m++) Zs(k, g, m) = static_cast<double>(draw[m]);
            }
        }
    }
}



// ---------------- Cosine similarity ----------------
inline double cosine_similarity(const arma::vec& x, const arma::vec& y) {
    double num = arma::dot(x, y);
    double denom = std::sqrt(arma::dot(x, x) * arma::dot(y, y));
    return (denom > 0) ? num / denom : 0.0;
}

// ---------------- update_P ----------------
arma::mat update_P(arma::mat P,
                   const arma::mat& Ap,
                   const arma::mat& Bp,
                   const std::vector<arma::cube>& Z_list,
                   const std::vector<arma::imat>& W_list,
                   const std::vector<arma::mat>& E_list,
                   const std::vector<arma::mat>& A_list,
                   int S,
                   const std::vector<int>& fix_colsP,
                   double psi = 1.0,
                   double similarity_threshold = 0.9,
                   int max_attempts = 10) {

    int i = P.n_rows;
    int n_signatures = P.n_cols;

    // 1️⃣ Atualiza P normalmente
    for (int m = 0; m < n_signatures; m++) {
        if (std::find(fix_colsP.begin(), fix_colsP.end(), m + 1) != fix_colsP.end()) continue;

        for (int k = 0; k < i; k++) {
            double shape = Ap(k, m) + 1.0;
            double rate  = Bp(k, m);

            for (int s = 0; s < S; s++) {
                const arma::cube& Zs = Z_list[s];
                const arma::imat& Ws = W_list[s];
                const arma::mat& Es  = E_list[s];
                const arma::mat& As  = A_list[s];

                shape +=   arma::accu(psi * As(m, m) * Zs.slice(m).row(k));
                rate  +=  psi * As(m, m) * arma::accu(Es.row(m)) * Ws(1, 1);
            }

            double scale = 1.0 /rate;           
            P(k, m) = R::rgamma(shape, scale);
        }
    }

    // 2️⃣ Verificação de similaridade e reamostragem
    for (int m1 = 0; m1 < n_signatures - 1; m1++) {
        for (int m2 = m1 + 1; m2 < n_signatures; m2++) {

            double curr_similarity = cosine_similarity(P.col(m1), P.col(m2));
            int attempt = 0;

            while (curr_similarity > similarity_threshold && attempt < max_attempts) {
                for (int k = 0; k < i; k++) {
                    double shape = Ap(k, m2) + 1.0;
                    double rate  = Bp(k, m2);

                    for (int s = 0; s < S; s++) {
                        const arma::cube& Zs = Z_list[s];
                        const arma::imat& Ws = W_list[s];
                        const arma::mat& Es  = E_list[s];
                        const arma::mat& As  = A_list[s];

                        shape +=  arma::accu(psi * As(m2, m2) * Zs.slice(m2).row(k));
                        rate  += psi * As(m2, m2) * arma::accu(Es.row(m2)) * Ws(1, 1);
                    }

                    double scale = 1.0 / rate;      
                    P(k, m2) = R::rgamma(shape, scale);
                }

                curr_similarity = cosine_similarity(P.col(m1), P.col(m2));
                attempt++;
            }
        }
    }

    return P;
}


// ----------------- update_P -----------------
arma::mat update_P_2(arma::mat P,
                   const arma::mat& Ap,
                   const arma::mat& Bp,
                   const std::vector<arma::cube>& Z_list,
                   const std::vector<arma::imat>& W_list,
                   const std::vector<arma::mat>& E_list,
                   const std::vector<arma::mat>& A_list,
                   int S,
                   const std::vector<int>& fix_colsP) {

  int i = P.n_rows;
  int n_signatures = P.n_cols;

  for (int m = 0; m < n_signatures; m++) {
    if (is_fixed(fix_colsP, m + 1)) continue;
    for (int k = 0; k < i; k++) {
      double shape = Ap(k, m) + 1.0;
      double rate  = Bp(k, m);
      for (int s = 0; s < S; s++) {
        const arma::cube& Zs = Z_list[s];// i x j x n_signatures
        const arma::imat& Ws = W_list[s];  
        const arma::mat& Es  = E_list[s]; // n x j
        const arma::mat& As  = A_list[s]; // n x n
        // CORREÇÃO: usar slice(m).row(k)
        shape +=  arma::accu(As(m, m) * Zs.slice(m).row(k));
        rate  += As(m, m) * arma::accu(Es.row(m))*Ws(1, 1);
      }
      double scale = 1.0 / rate;   
      P(k, m) = R::rgamma(shape, scale);
    }
  }
  return P;
}


// ----------------- update_E -----------------
std::vector<arma::mat> update_E(std::vector<arma::mat> E_list,
                                const std::vector<arma::mat>& Ae_list,
                                const std::vector<arma::mat>& Be_list,
                                const std::vector<arma::cube>& Z_list,
                                const std::vector<arma::imat>& W_list,
                                const arma::mat& P,
                                const std::vector<arma::mat>& A_list,
                                int S,
                                const std::vector<int>& fix_colsE) {

    int n_signatures = P.n_cols;

    for (int s = 0; s < S; s++) {

        arma::mat Es  = E_list[s];
        const arma::mat& Aes = Ae_list[s];
        const arma::mat& Bes = Be_list[s];
        const arma::cube& Zs = Z_list[s];
        const arma::imat& Ws = W_list[s];
        const arma::mat& As  = A_list[s];

        int J = Es.n_cols;

        for (int m1 = 0; m1 < n_signatures; m1++) {   // <- aqui muda

            for (int j = 0; j < J; j++) {

                if (is_fixed(fix_colsE, m1 + 1)) continue;

                double sumZ =  arma::accu(Zs.slice(m1).col(j));

                double sumP = As(m1, m1) * arma::accu(P.col(m1))*Ws(1, 1);

                double shape = Aes(m1, j) + 1.0 + As(m1, m1) * sumZ;

                double rate  = Bes(m1, j) + sumP;

                double scale = 1.0 / rate;

                Es(m1, j) = R::rgamma(shape, scale);
            }
        }

        E_list[s] = Es;
    }

    return E_list;
}


// ----------------- update_Bp_Be -----------------
void update_Bp_Be(arma::mat& Bp,
                  std::vector<arma::mat>& Be_list,
                  const arma::mat& P,
                  const std::vector<arma::mat>& E_list,
                  const arma::mat& Ap,
                  const std::vector<arma::mat>& Ae_list,
                  double ap, double ae, double bp, double be,
                  int S) {

    for (uword k = 0; k < P.n_rows; k++) {
        for (uword m = 0; m < P.n_cols; m++) {
            double shape = Ap(k, m) + 1.0 + ap;
            double rate  = P(k, m) + bp;
            double scale = 1.0 / rate;
            double val = R::rgamma(shape, scale);
            Bp(k, m) = val;
        }
    }

    for (int s = 0; s < S; s++) {
        arma::mat& Bes       = Be_list[s];
        const arma::mat& Es  = E_list[s];
        const arma::mat& Aes = Ae_list[s];
        for (uword m = 0; m < Es.n_rows; m++) {
            for (uword j = 0; j < Es.n_cols; j++) {
                double shape = Aes(m, j) + 1.0 + ae;
                double rate  = Es(m, j) + be;
                double scale =  1.0 / rate;
                double val = R::rgamma(shape, scale);
                Bes(m, j) = val;
            }
        }
    }
}


// ----------------- update_Ap_Ae -----------------
void update_Ap_Ae(arma::mat& Ap,
                  std::vector<arma::mat>& Ae_list,
                  const arma::mat& Bp,
                  const std::vector<arma::mat>& Be_list,
                  const arma::mat& P,
                  const std::vector<arma::mat>& E_list,
                  double lp,
                  double le,
                  int S) {

    // ===== Atualiza Ap =====
    for (uword k = 0; k < Ap.n_rows; k++) {
        // percorre as colunas em ordem aleatória (como sample no R)
        arma::uvec cols = arma::randperm(Ap.n_cols);
        for (uword idx = 0; idx < cols.n_elem; idx++) {
            uword m = cols[idx];

            double x = Ap(k, m);
            double y = R::rgamma(x, 1.0);

            double term_new = std::log(lp) - lp * y +
                ((y + 1.0) * std::log(Bp(k, m)) +
                 y * std::log(P(k, m)) - lgamma_R(y + 1.0));

            double term_old = std::log(lp) - lp * x +
                ((x + 1.0) * std::log(Bp(k, m)) +
                 x * std::log(P(k, m)) - lgamma_R(x + 1.0));

            double alpha = term_new - term_old +
                (R::dgamma(x, y, 1.0, true) - R::dgamma(y, x, 1.0, true));

            double rho = std::min(1.0, std::exp(alpha));
            if (!R_finite(rho)) rho = 0.0;

            if (R::runif(0.0, 1.0) <= rho)
                Ap(k, m) = y;
            else
                Ap(k, m) = x;
        }
    }

    // ===== Atualiza Ae =====
    for (int s = 0; s < S; s++) {
        arma::mat Aes = Ae_list[s];
        const arma::mat& Bes = Be_list[s];
        const arma::mat& Es  = E_list[s];

        arma::uvec rows = arma::randperm(Aes.n_rows);
        for (uword idx = 0; idx < rows.n_elem; idx++) {
            uword m = rows[idx];

            double x = Aes(m, 0);
            double y = R::rgamma(x, 1.0);

            // soma log(Be) e log(E) nas colunas (como sum(...) no R)
            double sum_new = 0.0, sum_old = 0.0;
            for (uword g = 0; g < Aes.n_cols; g++) {
                sum_new += ((y + 1.0) * std::log(Bes(m, g)) +
                            y * std::log(Es(m, g)) - lgamma_R(y + 1.0));
                sum_old += ((x + 1.0) * std::log(Bes(m, g)) +
                            x * std::log(Es(m, g)) - lgamma_R(x + 1.0));
            }

            double term_new = std::log(le) - le * y + sum_new;
            double term_old = std::log(le) - le * x + sum_old;

            double alpha = term_new - term_old +
                (R::dgamma(x, y, 1.0, true) - R::dgamma(y, x, 1.0, true));

            double rho = std::min(1.0, std::exp(alpha));
            if (!R_finite(rho)) rho = 0.0;

            if (R::runif(0.0, 1.0) <= rho)
                Aes.row(m).fill(y);
            else
                Aes.row(m).fill(x);
        }

        Ae_list[s] = Aes;
    }
}






// [[Rcpp::export]]
Rcpp::List burn_in_8(int num_iterations,
                     arma::mat P,
                     arma::mat Ap,
                     arma::mat Bp,
                     std::vector<arma::mat> E_list,
                     std::vector<arma::mat> Ae_list,
                     std::vector<arma::mat> Be_list,
                     std::vector<arma::mat> A_list,
                     std::vector<arma::cube> Z_list,
                     std::vector<arma::cube> Fi_list,
                     std::vector<arma::imat> M_list,
                     std::vector<arma::imat> W_list,
                     int S,
                     int i,
                     int n_signatures,
                     std::vector<int> j_list,
                     std::vector<int> fix_colsP,
                     std::vector<int> fix_colsA,
                     std::vector<int> fix_colsE,
                     double le,
                     double incremento_y,
                     double decressimo_l_1,
                     double decressimo_l_2,
                     double decressimo_l_3) {

    (void)i;

    // --- Valores fixos ---
     const double ap = 5.0, ae = 5.0, bp = 0.05, be = 0.01;   //<---------- Multi_Estudo

     // --- Valores fixos ---
    //const double ap = 0.8005, ae = 2.0506 , bp = 0.043, be = 1.325; //<---------- Hiper SigneR


    // --- Controle de psi (y) ---
    double y   = -2.0;
    double psi = std::pow(10.0, y);

    // --- Controle de l ---
    double l  = 1.0;
    double lp = 0.5;    //<---------- Multi_Estudo
    //double lp = 5.3329; //<---------- Hiper SigneR
    // --- Armazenamento ---
    int start_store = 7701;
    int end_store   = 13000;
    int n_store     = std::max(0, end_store - start_store + 1);

    arma::cube P_store(P.n_rows, P.n_cols, n_store, arma::fill::zeros);
    std::vector<std::vector<arma::mat>> A_store(S, std::vector<arma::mat>(n_store));
    std::vector<std::vector<arma::mat>> E_store(S, std::vector<arma::mat>(n_store));

    int store_idx = 0;

    // ======================
    // ==== LOOP PRINCIPAL ==
    // ======================
    for (int iter = 0; iter < num_iterations; iter++) {

        // =======================================================
        // 1) PSI — lógica idêntica ao R
        // =======================================================
        // R: if (iter %% 100 == 0)
        if ((iter + 1) % 100 == 0) {

            // R: if (y < 0) y <- y + incremento_y
            if (y < 0.0) {
                y += incremento_y;
            }

            // R: if (y >= 0) y <- 0
            if (y >= 0.0) {
                y = 0.0;
            }
        }

        // R: psi <- 10^y
        psi = std::pow(10.0, y);

        // =======================================================
        // 2) Ajuste de l (igual ao seu código)
        // =======================================================
        if ((iter + 1) % 101 == 0) {
            if (l > 0.20) {
                l -= decressimo_l_1;
            } else if (l > 0.01) {
                l -= decressimo_l_2;
            } else {
                l = std::max(l - decressimo_l_3, 0.0005);
            }
        }

        // =======================================================
        // 3) Atualizações principais
        // =======================================================
        if (iter < 7700) {
            A_list = update_A(P, E_list, A_list, M_list, W_list,
                              S, n_signatures, l, fix_colsA);
        } else {
            A_list = update_A_2(P, E_list, A_list, M_list, W_list,
                                S, n_signatures, l, fix_colsA);
        }

        update_Fi_Z(P, A_list, E_list, Z_list, Fi_list, M_list, S, n_signatures, j_list);

        if (iter < 7700) {
            P = update_P(P, Ap, Bp, Z_list, W_list, E_list, A_list, S, fix_colsP, psi);
        } else {
            P = update_P_2(P, Ap, Bp, Z_list, W_list, E_list, A_list, S, fix_colsP);
        }

        E_list = update_E(E_list, Ae_list, Be_list, Z_list, W_list, P, A_list, S, fix_colsE);

        update_Bp_Be(Bp, Be_list, P, E_list, Ap, Ae_list, ap, ae, bp, be, S);
       
        update_Ap_Ae(Ap, Ae_list, Bp, Be_list, P, E_list, lp, le, S);


        // =======================================================
        // 4) ARMAZENAMENTO
        // =======================================================
        if ((iter + 1) >= start_store && (iter + 1) <= end_store) {
            int save_pos = store_idx++;
            if (save_pos < n_store) {
                P_store.slice(save_pos) = P;
                for (int s = 0; s < S; s++) {
                    A_store[s][save_pos] = A_list[s];
                    E_store[s][save_pos] = E_list[s];
                }
            }
        }
    }

    // =======================================================
    // 5) RETORNO (SEM ALTERAÇÕES)
    // =======================================================
    return Rcpp::List::create(
        _["P"]          = P,
        _["E_list"]     = E_list,
        _["A_list"]     = A_list,
        _["Ap"]         = Ap,
        _["Ae_list"]    = Ae_list,
        _["Bp"]         = Bp,
        _["Be_list"]    = Be_list,
        _["Z_list"]     = Z_list,
        _["Fi_list"]    = Fi_list,
        _["M_list"]     = M_list,
        _["W_list"]     = W_list,
        _["P_store"]    = P_store,
        _["A_store"]    = A_store,
        _["E_store"]    = E_store,
        _["store_range"]= Rcpp::IntegerVector::create(start_store, end_store)
    );
}



// [[Rcpp::export]]
Rcpp::List burn_in_9(int num_iterations,
                     arma::mat P,
                     arma::mat Ap,
                     arma::mat Bp,
                     std::vector<arma::mat> E_list,
                     std::vector<arma::mat> Ae_list,
                     std::vector<arma::mat> Be_list,
                     std::vector<arma::mat> A_list,
                     std::vector<arma::cube> Z_list,
                     std::vector<arma::cube> Fi_list,
                     std::vector<arma::imat> M_list,
                     std::vector<arma::imat> W_list,
                     int S,
                     int i,
                     int n_signatures,
                     std::vector<int> j_list,
                     std::vector<int> fix_colsP,
                     std::vector<int> fix_colsA,
                     std::vector<int> fix_colsE,
                     double le,
                     double incremento_y,
                     double decressimo_l_1,
                     double decressimo_l_2,
                     double decressimo_l_3) {

    (void)i;
                     
    const double ap = 5.0, ae = 5.0, bp = 0.05, be = 0.01; //<---------- Hiper Multi-Estudo
     // --- Valores fixos ---
    //const double ap = 0.8005, ae = 2.0506 , bp = 0.043, be = 1.325; //<---------- Hiper SigneR
    double y = -2.0;
    double psi = std::pow(10.0, y);
    double l = 1.0;
    double lp = 0.5;  //<---------- Hiper Multi-Estudo
    //double lp = 5.3329; //<---------- Hiper SigneR
    // -------- Parâmetros de armazenamento --------
    int start_store = 1001;
    int end_store   = 5000;
    int n_store     = std::max(0, end_store - start_store + 1);

    arma::cube P_store(P.n_rows, P.n_cols, n_store, arma::fill::zeros);
    std::vector<std::vector<arma::mat>> A_store(S, std::vector<arma::mat>(n_store));
    std::vector<std::vector<arma::mat>> E_store(S, std::vector<arma::mat>(n_store));

    int store_idx = 0;

    // -------- Loop principal --------
    for (int iter = 0; iter < num_iterations; iter++) {

        // Ajuste do psi
        if ((iter + 1) % 100 == 0) {
            if (y < 0.0) y += incremento_y;
            if (y >= 0.0) y = 0.0;
        }
        psi = std::pow(10.0, y);

        // Ajuste de l
        if ((iter + 1) % 101 == 0) {
            if (l > 0.20)      l -= decressimo_l_1;
            else if (l > 0.01) l -= decressimo_l_2;
            else               l = std::max(l - decressimo_l_3, 0.0005);
        }

       
        A_list = update_A_2(P, E_list, A_list, M_list, W_list, S, n_signatures, l, fix_colsA);
       

        update_Fi_Z(P, A_list, E_list, Z_list, Fi_list, M_list, S, n_signatures, j_list);

        P = update_P_2(P, Ap, Bp, Z_list, W_list, E_list, A_list, S, fix_colsP);
      

        E_list = update_E(E_list, Ae_list, Be_list, Z_list, W_list, P, A_list, S, fix_colsE);

        update_Bp_Be(Bp, Be_list, P, E_list, Ap, Ae_list, ap, ae, bp, be, S);
        update_Ap_Ae(Ap, Ae_list, Bp, Be_list, P, E_list, lp, le, S);

        // ---------- Armazenamento ----------
        if ((iter + 1) >= start_store && (iter + 1) <= end_store) {
            int save_pos = store_idx++;
            if (save_pos < n_store) {
                P_store.slice(save_pos) = P;
                for (int s = 0; s < S; s++) {
                    A_store[s][save_pos] = A_list[s];
                    E_store[s][save_pos] = E_list[s];
                }
            }
        }
    }

    return Rcpp::List::create(
        _["P"]          = P,
        _["E_list"]     = E_list,
        _["A_list"]     = A_list,
        _["Ap"]         = Ap,
        _["Ae_list"]    = Ae_list,
        _["Bp"]         = Bp,
        _["Be_list"]    = Be_list,
        _["Z_list"]     = Z_list,
        _["Fi_list"]    = Fi_list,
        _["M_list"]     = M_list,
        _["W_list"]     = W_list,
        // ---- Traços armazenados ----
        _["P_store"]    = P_store,
        _["A_store"]    = A_store,
        _["E_store"]    = E_store,
        _["store_range"]= Rcpp::IntegerVector::create(start_store, end_store)
    );
}
