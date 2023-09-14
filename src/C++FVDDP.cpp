#include <Rcpp.h>
using namespace Rcpp;

//a function to comute the binomial coefficient
int binomial_int(int n, int k){

    //trivial case
    if (n==k){
        return 1;
    }

    //otherwise, try swapping to limit the size of the loop
    else{
        if (k < n-k){
            k = n-k;
        }

        //base for the result
        int res = 1;

        //looop through the values of k, on the numerator
        for(int i = 1; i<=n-k; i++){
            res *= (k+i);
        }

        //looop on the denominator
        for(int i = 1; i<=n-k; i++){
            res = res/i;
        }

        //return the value
        return res;
    }
}

//compute the probability C
double C_cpp(int m, int n, double t, NumericVector lambda, NumericMatrix& C_matrix){

    //if esult is not already stored in the matrix of results
    if(NumericVector::is_na(C_matrix(m,n))){

        //sex the sign of the result
        double r = pow(-1, m-n);

        //loop to multiply the left hand side product
        for(int j = n+1; j <= m; j++){
            r *= lambda[j];
        }

        //initilize the right hand side sum
        double s = 0;

        //it is a sum of products
        for(int k = n; k <= m; k++){

            //make it start from one
            double p = 1;

            //loop to get the denominator
            for(int j = n; j <= m; j++){

                //multiplying each time
                if (j != k){
                    p *= lambda[k] - lambda[j];
                }
            }

            //recall the exponential term on the numerator
            s += exp(-t*lambda[k])/p;
        }

        //get the final result and save it in the matrix and return it
        double res = r*s;
        C_matrix(m,n) = res;
        return res;
    }

    //otherwise, just return the value
    else{
        return C_matrix(m,n);
    }
}

//compute the transsition probability p
double p_cpp(NumericVector v_m, NumericVector v_n, double t, NumericVector lambda, NumericMatrix &C_matrix){

    //iniilaize some values
    int m = sum(v_m);
    int n = sum(v_n);
    double p_t;

    //in case the vectos are equal
    if(n == m){
        p_t = exp(-t * lambda[m]);
    }

    //if not
    else{

        //comute the numeratorof the densit of a multivariate hypergometric
        NumericVector vect_binom = choose(v_m, v_n);
        int prod = std::accumulate(vect_binom.begin(), vect_binom.end(), 1, std::multiplies<int>());

        //merge it with C and the denominator to get p
        p_t = C_cpp(m, n, t, lambda, C_matrix) * prod / binomial_int(m, n);
    }

    //the result
    return p_t;
}


// [[Rcpp::export]]
NumericVector compute_new_weights_cpp(NumericMatrix M, NumericMatrix LM, double t, NumericVector w, NumericVector lambda, NumericMatrix& C_matrix){

    //the emptu new vector of weights
    NumericVector w_new(LM.nrow());

    //loop though it
    for(int i = 0; i < LM.nrow(); i++){

        //print the stte of the procedure
        Rcout << "iteration " << i+1 << " of " << LM.nrow();

        //he weight is 0 by defaul
        double weight = 0;

        //take the row of L(M) as a vector
        NumericVector v_LM = LM(i,_);

        //this boolean variable tells if there are ancetors in M; by default, no
        bool is_in_LM = false;

        //iterate throught the rows of M
        for (int q = 0; q < M.nrow(); q++){

            //get it as a vector
            NumericVector v_M = M(q,_);

            //check whether v_LM <= v_M
            bool is_leq_v_M = true;

            //iterate on both the vectors
            for(int j = 0; j < M.ncol(); j++){

                //if the condition is broken, the brevious boolean is false
                if(v_LM[j] > v_M[j]){
                    is_leq_v_M = false;
                    break;
                }
            }

            //if v_M is an ancestor of v_LM
            if(is_leq_v_M == true){

                //the vector v_LM is in L(M)
                if(is_in_LM == false){
                    is_in_LM = true;
                }

                //compute the transition probability
                double p_t = p_cpp(v_M, v_LM, t, lambda, C_matrix);

                //add it to the weight
                weight += p_t * w[q];
            }
        }

        //hence if the vector is in LM, save its weight
        if(is_in_LM == true){
            w_new[i] = weight;
        }

        //otherwise save a NA, to be revode later
        else{
            w_new[i] = NA_REAL;
        }

        //clear console line
        Rcout << "\r";
    }

    //new console line
    Rcout << std::endl;

    //return the vector with weights
    return w_new;
}


//compute the cartesian product of a list of vector, as a matrix
NumericMatrix cartesian_cpp(List l){

    //initialize some values
    int N = 1;
    int K = l.length();

    //N is the number of rows
    for(int j = 0; j < K; j++){
        NumericVector v = l[j];
        N *= v.length();
    }

    //crate the empty matrix
    NumericMatrix M(N, K);

    //numer of replicates of each value in row vectors
    int reps = 1;

    //loop to mimic all possible combinations, from last row
    for(int j = K-1; j >= 0; j--){

        //take the vector
        NumericVector v = l[j];

        //change the number of replicates
        reps *= v.length();

        //repeat it and save it as a column of the matrix
        M(_, j) = rep(rep_each(v, N/reps), reps);
    }

    //return the matrix
    return M;
}

//compute the marginal int the atomic case
double marginal_cpp(NumericVector v, double theta, NumericVector theta_P0_ystar,
std::map<std::vector<int>, double> &m_map){

    //save the vector as a key
    std::vector<int> key = as<std::vector<int>>(v);

    //check if it already is in te map
    if(m_map.count(key) == 1){
        return m_map[key];
    }

    //otherwise it has to be computed
    else{

        //the firs term, a Pochhammer's symbol
        double m = std::tgamma(theta)/std::tgamma(theta + sum(v));

        //loop on the elements of the vectors
        for(int i = 0; i < v.length(); i++){

            //is there are positive multipliciies, compute the new Pochhammer's term
            if(v[i] != 0){
                m *= std::tgamma(theta_P0_ystar[i] + v[i])/std::tgamma(theta_P0_ystar[i]);
            }
        }

        //sae the value in the map
        m_map[key] = m;

        //return it
        return m;
    }
}

//the shared component in the nonatomic case
double shared_w(int sum_k_past, int sum_n, int sum_k_future, double theta,
NumericMatrix &nonatomic_w_matrix){

    //if the value has not beeen computed
    if(NumericVector::is_na(nonatomic_w_matrix(sum_k_past, sum_k_future))){

        //cpompute both the argumets of the gamma for numerator and for the ndenominator
        NumericVector num = {theta + sum_k_past, theta + sum_n, theta + sum_k_future};
        NumericVector den = {theta + sum_k_past + sum_n + sum_k_future, theta, theta};

        //actually pass through th gamma and divide elementwise
        NumericVector res = gamma(num)/gamma(den);

        //get the product
        double u = std::accumulate(res.begin(), res.end(), 1.0, std::multiplies<double>());

        //save it in the appropriate place
        nonatomic_w_matrix(sum_k_past, sum_k_future) = u;

        //return the value
        return u;
    }

    //otherwise take it and return it
    else{
        return nonatomic_w_matrix(sum_k_past, sum_k_future);
    }

}

//generate the set D as a matrix, and compute its vector of weights
List generate_D_cpp(NumericVector n_past, NumericVector n_future, NumericVector n, double t_past,
double t_future, bool atomic, NumericVector lambda, NumericMatrix &C_matrix_past,
NumericMatrix &C_matrix_future, double theta, NumericVector theta_P0_ystar,
std::map<std::vector<int>, double> &m_map, NumericMatrix &nonatomic_w_matrix){

    //specify the vector size
    int K = n.length();

    //this will be necessay to identify shaed components
    NumericVector base(2*K);

    //this vector stores the indices of hared components
    std::vector<int> S_vect;

    //if the measure is nonatomic, take into account shared types among times
    if(atomic == false){

        //loop on the vectors
        for(int i=0; i < K; i++){

            //shared values from the past; in case start by 1
            if((n_past[i] >0) && ((n[i] > 0) || (n_future[i] > 0))){
                    base[i] = 1;
            }

            //shared values from the future; in case start by 1
            if((n_future[i] >0) && ((n[i] > 0) || (n_past[i] > 0))){
                    base[i+K] = 1;
            }

            //if there is sharing, append the index to S
            if((base[i] == 1) || (base[i+K])){
                S_vect.push_back(i);
            }
        }
    }

    //change its type
    NumericVector S = wrap(S_vect);

    //initialize an empty list
    List l(2*K);

    //loop to create the vectors
    for(int i=0; i < K; i++){

        //thi is for the past, save it in the list
        IntegerVector v_past = seq(base[i], n_past[i]);
        l[i] = v_past;

        //this is for the future, save it in the list too
        IntegerVector v_future = seq(base[i+K], n_future[i]);
        l[i+K] = v_future;
    }

    //get the mtrix using the cartesian product
    NumericMatrix D = cartesian_cpp(l);

    //initialize an empty vector
    NumericVector w(D.nrow());

    //if there is just one row, the weight is 1
    if(D.nrow() == 1){
        w={1};
    }

    //othewise
    else{

        //loop on the rows, each is a pair of vectors
        for(int i = 0; i < D.nrow(); i++){

            //take a line and get the the two vectors
            NumericVector line = D(i, _);
            NumericVector k_past = line[Range(0,K-1)];
            NumericVector k_future = line[Range(K, 2*K-1)];

            //compute the transition probabilities
            double p_past = p_cpp(n_past, k_past, t_past, lambda, C_matrix_past);
            double p_future = p_cpp(n_future, k_future, t_future, lambda, C_matrix_future);

           //initialize the weights
            double u;

            //in the atomic case
            if(atomic == true){

                //compute all terms of the formula
                double m_sum = marginal_cpp(k_past + n + k_future, theta, theta_P0_ystar, m_map);
                double m_past = marginal_cpp(k_past, theta, theta_P0_ystar, m_map);
                double m_n = marginal_cpp(n, theta, theta_P0_ystar, m_map);
                double m_future = marginal_cpp(k_future, theta, theta_P0_ystar, m_map);

                //merge them together to get the expression
                u = m_sum/(m_past * m_n * m_future);
            }

            //in the nonatomic case
            if(atomic == false){

                //save its first component
                u = shared_w(sum(k_past), sum(n), sum(k_future), theta, nonatomic_w_matrix);

                //if thre are shared types, compute the right hand side of the formula
                if(S.length() != 0){

                    //to do it, iterate through all indices
                    for(int s = 0; s < S.length(); s++){

                        //compute each term
                        double gamma_sum = std::tgamma(k_past[S[s]] + n[S[s]] + k_future[S[s]]);
                        double gamma_past = (k_past[S[s]] > 0) ? std::tgamma(k_past[S[s]]) : 1;
                        double gamma_n = (n[S[s]] > 0) ? std::tgamma(n[S[s]]) : 1;
                        double gamma_future = (k_future[S[s]] > 0) ? std::tgamma(k_future[S[s]]) : 1;

                        //merge them together, to get the result
                        u *= gamma_sum/(gamma_past * gamma_n *gamma_future);
                    }
                }
            }

            //save the unnormalized weight
            w[i] = p_past * p_future * u;
        }

       //normalize the vector of weights
        w = w/sum(w);
    }

    //return the list
    return List::create(_["D"] = D, _["w"] = w);
}

// [[Rcpp::export]]
List compute_M_w_cpp(NumericMatrix M_past, NumericMatrix M_future, NumericVector n, double t_past, double t_future,
NumericVector w_past, NumericVector w_future, bool atomic, NumericVector lambda,
NumericMatrix &C_matrix_past, NumericMatrix &C_matrix_future, double theta, NumericVector theta_P0_ystar,
NumericMatrix &nonatomic_w_matrix){

    //this maps takes a vector pair as a key, and stores their weight
    std::map<std::vector<int>, double> weights;

    //this map takes a vector as key, and saves the value of the nonatomic marginal
    std::map<std::vector<int>, double> m_map;

    //save the value of K
    int K = M_past.ncol();

    //iterate on past multiplicities vector
    for (int i_past = 0; i_past < M_past.nrow(); i_past++){
        NumericVector n_past = M_past(i_past, _);

        //iterate on future multiplicities vectors
        for(int i_future = 0; i_future < M_future.nrow(); i_future++){
            NumericVector n_future = M_future(i_future, _);

            //prin thte state of the system
            Rcout << "iteration " << i_past * M_future.nrow() + i_future + 1 << " out of "
            << M_past.nrow() * M_future.nrow();

            //generate the set D, with its weights
            List D_w = generate_D_cpp(n_past, n_future, n, t_past, t_future, atomic, lambda, C_matrix_past,
            C_matrix_future, theta, theta_P0_ystar, m_map, nonatomic_w_matrix);

            //store the elements of the list
            NumericMatrix D = D_w[0];
            NumericVector w = D_w[1];

            //loop on the rows of D
            for(int i = 0; i < D.nrow(); i++){

                //save a line and split it itno two vectors, for past and future
                NumericVector line = D(i, _);
                NumericVector k_past = line[Range(0,K-1)];
                NumericVector k_future = line[Range(K, 2*K-1)];

                //the result to append is their sum, with n too
                NumericVector key = k_past + k_future + n;

                //update adding the weight to the approprite vector
                weights[as<std::vector<int> >(key)] += w[i] * w_past[i_past] * w_future[i_future];
            }

            //new conole line
            Rcout << "\r";
        }
    }

    //print the state
    Rcout << std::endl << "Merging and sorting in progress... ";

    //finally get the mtrix and the vector
    NumericVector w_new (weights.size());
    NumericMatrix M_new (weights.size(), K);

    //create an  iterator on the map
    int j = 0;
    std::map<std::vector<int>, double>::iterator it;

    //iterate on the map
    for(it = weights.begin(); it!=weights.end(); ++it){

        //take the key and store it as a row of the matrix
        NumericVector v = wrap(it->first);
        M_new(j,_) = v;

        //take the value and save it in the vector of weights
        w_new[j] = it->second;
        j++;
    }

    //print that the computation is over
    Rcout << "Completed!" << std::endl;

    //return both the matrix and its weight
    return List::create(_["M"] = M_new , _["w"] = w_new);
}

//generate a vector via a multivariate hypergeometric distribution
NumericVector rMVH_cpp(NumericVector n, int k){

    //trivial case
    if(k == sum(n)){
        return n;
    }

    //the same
    if(k == 0){
        return NumericVector (n.length());
    }

    //progressively remove values
    else{

        //initialize the output, and clone the other vector not to modify it
        NumericVector m(n.length());
        NumericVector p = clone(n);

        //repeat k times
        for(int i= 1; i<=k; i++){

            //sample a new value
            int j = sample(n.length(), 1, false, p)[0];

            //shif it from a vector to the other
            p[j-1] -= 1;
            m[j-1] += 1;
        }

        //return the vector built this way
        return m;
    }
}


//simulate the decay via a 1-dimensional process
NumericVector decay_cpp(NumericVector v, double t, NumericVector lambda){

    //set the starting cardinality
    int m = sum(v);

    //repeatedly draw wxponential and reduce time
    while(m > 0){
        t -= rexp(1, lambda[m])[0];

        //if time is over, stop
        if (t <= 0){
            break;
        }

        //otherwise reduce the cardinality
        m -= 1;
    }

    //get a multidimensional vector drawing from multivariate hypergeomentric
    NumericVector v_LM = rMVH_cpp(v, m);

    //return the vector
    return v_LM;
}

// [[Rcpp::export]]
List montecarlo_sample_prop_cpp(NumericMatrix M, double t, int N, NumericVector w, NumericVector lambda){

    //create a map of weights
    std::map<std::vector<int>, int> weights;

    //repeat this for each sample drawn
    for(int j = 1; j <= N; j++){

        //print the percentage
        Rcout << std::fixed << std::setprecision(2) << j*100/N << "%  ";

        //select a row of M with probability given by w
        int i = sample(M.nrow(), 1, false, w)[0] - 1;

        //save it as a vector
        NumericVector v_M = M(i,_);

        //make it decay in the given time
        NumericVector v_LM = decay_cpp(v_M, t, lambda);

        //increase by one the weight of such vector
        weights[as<std::vector<int>>(v_LM)] += 1;

        //clear console line
        Rcout << "\r";

    }

    //initialize the output, and clone the other vector not to modify it
    NumericVector w_new (weights.size());
    NumericMatrix LM (weights.size(), M.ncol());

    //create an iterator for a map, and an iteger to specify indices
    int j = 0;
    std::map<std::vector<int>, int>::iterator it;

    //iterate through the map
    for(it = weights.begin(); it!=weights.end(); ++it){

        //take the key and save it as a row of the matrix
        NumericVector v_LM = wrap(it->first);
        LM(j,_) = v_LM;

        //take the value and save it in the vector of weights
        w_new[j] = it->second;
        j++;
    }

    //return both the matrixx and its unnormalized weight
    return List::create(_["LM"] = LM , _["w"] = w_new);
}


//compute the importance of a pair, given the ancestors
double importance_cpp(NumericVector k_past, NumericVector k_future, NumericVector n, NumericVector n_past,
NumericVector n_future, bool atomic, double theta, NumericVector theta_P0_ystar,
std::map<std::vector<int>, double> &m_map, NumericMatrix &nonatomic_w_matrix){

    //set the length of vectors
    int K = n.length();

    //initialize the importance
    double imp;

    //in case the measure is atomic
    if(atomic == true){

        //compute the marginal for each vector
        double m_sum = marginal_cpp(k_past + n + k_future, theta, theta_P0_ystar, m_map);
        double m_past = marginal_cpp(k_past, theta, theta_P0_ystar, m_map);
        double m_n = marginal_cpp(n, theta, theta_P0_ystar, m_map);
        double m_future = marginal_cpp(k_future, theta, theta_P0_ystar, m_map);

        //use the formula to compute the imporance
        imp = m_sum/(m_past * m_n * m_future);
    }

    //otherwise
    else if(atomic == false){

        //initialize the vector of shared values
        std::vector<int> S_vect;

        //iterate through the vectors
        for(int i = 0; i < K; i++){

            //this boolean tells if there is information sharing for such index
            bool push = false;

            //check from the past
            if(((n_past[i] >0) && ((n[i] > 0) || (n_future[i] > 0)))  ){

                //if the vector is not in D return 0
                if(k_past[i] == 0){
                    return 0;
                }

                //note that  are shared types
                push = true;
            }

            //check from the future
            if(((n_future[i] >0) && ((n[i] > 0) || (n_past[i] > 0)))  ){

                //if the vector is not in D return 0
                if(k_future[i] == 0){
                    return 0;
                }

                //note that there are shared types
                push = true;
            }

            //if there are shared type, append the index to S
            if(push == true){
                S_vect.push_back(i);
            }
        }

        //convert S
        NumericVector S = wrap(S_vect);

        //compute the left hand side of the formula for the nonatomic case
        imp = shared_w(sum(k_past), sum(n), sum(k_future), theta, nonatomic_w_matrix);

        //if there is information sharing among times
        if(S.length() != 0){

            //loop on the indices
            for(int s = 0; s < S.length(); s++){

                //compute the right hand side, made by factorials
                double gamma_sum = std::tgamma(k_past[S[s]] + n[S[s]] + k_future[S[s]]);
                double gamma_past = (k_past[S[s]] > 0) ? std::tgamma(k_past[S[s]]) : 1;
                double gamma_n = (n[S[s]] > 0) ? std::tgamma(n[S[s]]) : 1;
                double gamma_future = (k_future[S[s]] > 0) ? std::tgamma(k_future[S[s]]) : 1;

                //mutliply it times the right hand side to get the importance
                imp *= gamma_sum/(gamma_past * gamma_n *gamma_future);
            }
        }
    }

    //return the value computed
    return imp;
}

//[[Rcpp::export]]
List montecarlo_sample_smooth_cpp(NumericMatrix M_past, NumericMatrix M_future, NumericVector n, double t_past,
double t_future, NumericVector w_past, NumericVector w_future, int N, bool atomic, NumericVector lambda,
double theta, NumericVector theta_P0_ystar, NumericMatrix &nonatomic_w_matrix){

    //se the length of vectors
    int K = n.length();

    //initialize a map for the weights, whose keys are vector pairs
    std::map<std::vector<int>, double> weights;

    //initialize a map of bins, i.e. from pair of ancestors vector to pairs of offsping vectors
    std::map<std::vector<int>, std::map<std::vector<int>, NumericVector>> bins;

    //initialize a map fr the marginal
    std::map<std::vector<int>, double> m_map;

    //repeat the amount of times of samples to draw
    for(int j = 1; j <= N; j++){

        //print the percentage
        Rcout << std::fixed << std::setprecision(2) << j*100/N << "%  ";

        //get a row from the past and from the future according tho the weights
        int i_past = sample(M_past.nrow(), 1, false, w_past)[0] - 1;
        int i_future = sample(M_future.nrow(), 1, false, w_future)[0] - 1;

        //save the ancestors vectors
        NumericVector n_past = M_past(i_past,_);
        NumericVector n_future = M_future(i_future,_);

        //simulate the decay to get offspring vectors
        NumericVector k_past = decay_cpp(n_past, t_past, lambda);
        NumericVector k_future = decay_cpp(n_future, t_future, lambda);

        //save the ancestor together to compose a key for their bin
        NumericVector bin_key (2*K);
        bin_key[Range(0, K-1)] = n_past;
        bin_key[Range(K, 2*K-1)] = n_future;

        //save the ancestors together to get an item of the bin
        NumericVector k_pair (2*K);
        k_pair[Range(0, K-1)] = k_past;
        k_pair[Range(K, 2*K-1)] = k_future;

        //if the bin does not exist or the item is not in the bin
        if ((bins.count(as<std::vector<int>>(bin_key)) == 0) || (bins[as<std::vector<int>>(bin_key)].count(as<std::vector<int>>(k_pair)) == 0)){

            //compute the importance of the offsping pair given the ancestors
            double imp = importance_cpp(k_past, k_future, n, n_past, n_future, atomic, theta, theta_P0_ystar, m_map, nonatomic_w_matrix);

            //include the count and the importance in the bin, mapped by the offsping pair
            bins[as<std::vector<int>>(bin_key)][as<std::vector<int>>(k_pair)] = NumericVector {imp, 1};
        }

        //ortherwise, is the bin and the pair exist
        else{

            //enter the bin, find the key offsping pair, increase its count by 1
            bins[as<std::vector<int>>(bin_key)][as<std::vector<int>>(k_pair)][1] += 1;
        }

        //clear the line
        Rcout << "\r";
    }

    //create an itherstor from a map to a map
    std::map<std::vector<int>, std::map<std::vector<int>, NumericVector>>::iterator bin_it;

    //use it to iterate through the bins
    for(bin_it = bins.begin(); bin_it != bins.end(); ++bin_it){

        //initialize the sum of the vector
        double w_sum = 0;

        //save the bin as a map
        std::map<std::vector<int>, NumericVector> bin = bin_it->second;

        //create an itherator from a vector pair to a numeric vector
        std::map<std::vector<int>, NumericVector>::iterator it_pair;

        //use it to iterate on each bin
        for(it_pair = bin.begin(); it_pair != bin.end(); ++it_pair){

           //in paricular, compute tje sum of importances
            w_sum += (it_pair->second)[0];
        }

        //if there is no importance on the bin, move on
        if(w_sum == 0){
            continue;
        }

        //iterate once more within the bin
        for(it_pair = bin.begin(); it_pair != bin.end(); ++it_pair){

            //then get the values of such pair, importance and count
            NumericVector val = it_pair->second;

            //if the pair has null importance, do not add it
            if(val[0] == 0){
                continue;
            }

            //get the vector pair
            NumericVector v = wrap(it_pair->first);

            //use as key for the weight the sum of the pair and n
            NumericVector key = v[Range(0, K-1)] + n + v[Range(K, 2*K-1)];

            //use them to update the weight, normalizing the constant
            weights[as<std::vector<int>>(key)] += (val[0]/w_sum) * val[1];
        }
    }

    //initialize an empty vector and matrix
    NumericVector w_new (weights.size());
    NumericMatrix M_new (weights.size(), K);

    //initialize an iterator and and a counter variable
    int j = 0;
    std::map<std::vector<int>, double>::iterator it;

    //use the iterrator on the map of weights
    for(it = weights.begin(); it!=weights.end(); ++it){

        //get the key vector and append it as a row
        NumericVector v_M = wrap(it->first);
        M_new(j,_) = v_M;

        //get the value and save it as the weight
        w_new[j] = it->second;
        j++;
    }

    //return a list with both the matrix and the weight
    return List::create(_["M"] = M_new , _["w"] = w_new);
}
