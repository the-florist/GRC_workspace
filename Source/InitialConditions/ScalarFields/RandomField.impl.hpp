/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(RANDOMFIELD_HPP_)
 #error "This file should only be included via RandomField.hpp"
 #endif

 #ifndef RANDOMFIELD_IMPL_HPP_
 #define RANDOMFIELD_IMPL_HPP_

 inline RandomField::RandomField(params_t a_params, InitialScalarData::params_t a_bkgd_params, std::string a_spec_type)
    : m_params(a_params), m_bkgd_params(a_bkgd_params), m_spec_type(a_spec_type)
{
    const double Mp = 1./m_bkgd_params.E;
    kstar = M_PI*((double) m_params.Nf)/m_params.L; //50.*(2.*M_PI/m_params.L);
    epsilon = m_params.L/30.; //100.;//0.5;//0.25 * (sqrt(3.)*2.*M_PI/m_params.L); //0.5;
    H0 = sqrt((8.0 * M_PI/3.0/pow(Mp, 2.))
            * (0.5*m_bkgd_params.velocity*m_bkgd_params.velocity 
                + 0.5*pow(m_bkgd_params.m * m_bkgd_params.amplitude, 2.0)));
    
    norm = m_params.A * pow(2.*M_PI/m_params.L, 3.);

    calc_spectrum();
}


template <class data_t>
void RandomField::compute(Cell<data_t> current_cell) const
{
    // Pull out the grid parameters
    int Nc = m_params.N;
    int N = m_params.Nf;
    int skip = (int)(N/Nc);

    double L = m_params.L;
    double dx = L/Nc;

    // Setting the lut that maps polarisation vectors to 
    // polarisation tensors.
    int lut[3][3];  
    lut[0][0] = 0;
    lut[0][1] = 1;
    lut[0][2] = 2;
    lut[1][0] = 1;
    lut[1][1] = 3;
    lut[1][2] = 4;
    lut[2][0] = 2;
    lut[2][1] = 4;
    lut[2][2] = 5;

    Coordinates<data_t> coords(current_cell, dx, m_params.center);

    // Coordinate of this cell in program units
    data_t x = coords.x + L/2;
    double y = coords.y + L/2;
    double z = coords.z + L/2;

    // Coordinates of this cell, unitless
    int i = static_cast<int>(x / dx);
    int j = static_cast<int>(y / dx);
    int k = static_cast<int>(z / dx);

    // This is to guard against ghost cells that can take you outside 
    // the domain of dependence of the box. Uses periodic BCs.
    if(i < 0) { i = Nc + i; }
    else if(i >= Nc) { i = i - Nc; }

    if(j < 0) { j = Nc + j; }
    else if(j >= Nc) { j = j - Nc; }

    if(k < 0) { k = Nc + k; }
    else if(k >= Nc) { k = k - Nc; }

    int rc = k + Nc*(j + Nc*i);

    // Error trap, to make sure we stay in the domain of dependence
    if (rc < 0)
    {
        cout << rc << endl;
        MayDay::Error("RandomField: Cell index value below zero.");
    }
    else if(rc > pow(Nc, 3.))
    {
        cout << rc << endl;
        MayDay::Error("RandomField: Cell index greater than Nc^3 at coarsest level.");
    }

    // The flattened position (leading with z or k)
    int r = (k + N*(j + N*i))*skip;

    // Error trap, to make sure we stay in the domain of dependence of the larger grid
    if (r > pow(N, 3.))
    {
        cout << r << endl;
        MayDay::Error("RandomField: Cell index greater than N^3 at coarsest level.");
    }

    // Assign position and momenum to the vars. objects
    double trace = 0;
    for(int l=0; l<3; l++) for(int p=l; p<3; p++) 
    {
        hx[lut[l][p]][r] *= norm; 

        if(m_spec_type == "position")
        {
            if (l==p) { hx[lut[l][p]][r] += 1.; }
            trace = abs(hx[0][r] + hx[3][r] + hx[5][r] - 3.);
        }
        else if(m_spec_type == "velocity")
        {
            hx[lut[l][p]][r] *= -1.;
            trace = abs(hx[0][r] + hx[3][r] + hx[5][r]);
        }
        else { MayDay::Error("Spectral type provided is an invalid option."); }
    } 

    // Trace free test
    if(trace > 1.e-12) 
    { 
        std::cout << fixed << setprecision(12);
        std::cout << "Field: " << m_spec_type << "\n";
        std::cout << "Trace of field is large here: \n";
        std::cout << "(" << i << "," << j << "," << k << ")\n";
        std::cout << rc << "," << r << "\n";
        std::cout << hx[0][r] << "," << hx[3][r] << "," << hx[5][r] << "\n";
        std::cout << trace << "\n";
    }

    if(m_spec_type == "position")
    {
        current_cell.store_vars(hx[lut[0][0]][r], c_h11);
        current_cell.store_vars(hx[lut[0][1]][r], c_h12);
        current_cell.store_vars(hx[lut[0][2]][r], c_h13);
        current_cell.store_vars(hx[lut[1][1]][r], c_h22);
        current_cell.store_vars(hx[lut[1][2]][r], c_h23);
        current_cell.store_vars(hx[lut[2][2]][r], c_h33);
    }
    else if(m_spec_type == "velocity")
    {
        current_cell.store_vars(hx[lut[0][0]][r], c_A11);
        current_cell.store_vars(hx[lut[0][1]][r], c_A12);
        current_cell.store_vars(hx[lut[0][2]][r], c_A13);
        current_cell.store_vars(hx[lut[1][1]][r], c_A22);
        current_cell.store_vars(hx[lut[1][2]][r], c_A23);
        current_cell.store_vars(hx[lut[2][2]][r], c_A33);
    }
    else { MayDay::Error("Spec type provided is invalid."); }
}

inline void RandomField::clear_data()
{
    pout() << "Clearing memory allocated to hx array.\n";
    for(int s=0; s<6; s++) { free(hx[s]); } // This causes a seg fault as is
    free(hx);
}

inline void RandomField::calc_spectrum()
{
    int N = m_params.Nf;
    int which_seed = 3;

    // Setting the lut that maps polarisation vectors to 
    // polarisation tensors.
    int lut[3][3];  
    lut[0][0] = 0;
    lut[0][1] = 1;
    lut[0][2] = 2;
    lut[1][0] = 1;
    lut[1][1] = 3;
    lut[1][2] = 4;
    lut[2][0] = 2;
    lut[2][1] = 4;
    lut[2][2] = 5;

    // Polarisation basis vectors and k
    double mhat[3] = {0., 0., 0.};
    double nhat[3] = {0., 0., 0.}; 
    double kmag = 0.;

    // Allocate memory for hij, mode functions and plans
    fftw_complex** hk;
    hk = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*) * 6);
    fftw_plan hij_plan[6];

    hx = (double**) malloc(sizeof(double*) * 6);
    for(int l=0; l<6; l++)
    {
        hk[l] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
        hx[l] = (double*) malloc(sizeof(double) * N * N * N);
        hij_plan[l] = fftw_plan_dft_c2r_3d(N, N, N, hk[l], hx[l], FFTW_ESTIMATE);
    }

    fftw_complex (*hplus);
    fftw_complex (*hcross);
    hplus = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
    hcross = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));

    // Extra memory for reality check on h+ (only for debug)
    double (*hplusx);
    hplusx = (double*) malloc(sizeof(double) * N * N * N);
    double (*hcrossx);
    hcrossx = (double*) malloc(sizeof(double) * N * N * N);

    fftw_plan plan1 = fftw_plan_dft_c2r_3d(N, N, N, hplus, hplusx, FFTW_ESTIMATE);
    fftw_plan plan2 = fftw_plan_dft_c2r_3d(N, N, N, hcross, hcrossx, FFTW_ESTIMATE);
    
    // Set all arrays to 0
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        if(k <= N/2)
        {
            for(int s=0; s<2; s++)
            {
                hplus[k + (N/2+1)*(j + N*i)][s] = 0.;
                hcross[k + (N/2+1)*(j + N*i)][s] = 0.;
            }
        }
        hplusx[k + N*(j + N*i)] = 0.;
        hcrossx[k + N*(j + N*i)] = 0.;

        for (int m=0; m<6; m++)
        {
            if(k <= N/2) { for(int s=0; s<2; s++) { hk[m][k + (N/2+1)*(j + N*i)][s] = 0.; } }
            hx[m][k + N*(j + N*i)] = 0.;
        }
    }

    // Set up random number generators (one independent seed per random draw)
    std::vector<int> seeds(10, 0);
    seeds[0] = 3539263;
    seeds[1] = 7586572;
    seeds[2] = 5060982;
    seeds[3] = 6793957;
    seeds[4] = 4764135;
    seeds[5] = 2034988;
    seeds[6] = 9635753;
    seeds[7] = 9350886;
    seeds[8] = 6855322;
    seeds[9] = 2933414;

    int seed = seeds[which_seed];
    default_random_engine engine(seed);
    uniform_real_distribution<double> theta_dist(0, 2*M_PI);
    uniform_real_distribution<double> sigma_dist(0, 1);
    double plus_mod, cross_mod, plus_arg, cross_arg;

    pout() << "Starting RandomField loop for " << m_spec_type << " field.\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
    {
        // Putting the 0 mode in the right spot in memory
        int I = i;
        int J = j;
        if(i > N/2) { I = invert_index(i, N); }
        if(j > N/2) { J = invert_index(j, N); }

        kmag = (double)(pow((pow((double)I, 2.0) + pow((double)J, 2.0) + pow((double)k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5));

        plus_mod = sigma_dist(engine);
        cross_mod = sigma_dist(engine);
        plus_arg = theta_dist(engine);
        cross_arg = theta_dist(engine);

        // Start of with random numbers filling the entire array
        // Real parts of h+, hx and hij
        if(kmag != 0)
        {
            if(i==4 && j==0 && k==0)
            {
                if(m_spec_type == "position") { 
                    hplus[k + (N/2+1)*(j + N*i)][0] = 0.;//1.;
                    hcross[k + (N/2+1)*(j + N*i)][0] = 0.;//1.;

                    /*hplus[k + (N/2+1)*(j + N*i)][0] = 1./sqrt(2.);
                    hcross[k + (N/2+1)*(j + N*i)][0] = 1./sqrt(2.);

                    hplus[k + (N/2+1)*(j + N*i)][1] = 1./sqrt(2.);
                    hcross[k + (N/2+1)*(j + N*i)][1] = 1./sqrt(2.);*/
                }
                else if(m_spec_type == "velocity")
                {
                    hplus[k + (N/2+1)*(j + N*i)][0] = 1.;//kmag;
                    hcross[k + (N/2+1)*(j + N*i)][0] = 1.;//kmag;
                }
            }

            for(int s=0; s<2; s++)
            {
                //hplus[k + (N/2+1)*(j + N*i)][s] = sqrt(-2. * log(plus_mod) * find_rayleigh_factor(kmag, m_spec_type));
                //hcross[k + (N/2+1)*(j + N*i)][s] = sqrt(-2. * log(cross_mod) * find_rayleigh_factor(kmag, m_spec_type));

                //if(s==0) { hplus[k + (N/2+1)*(j + N*i)][s] *= cos(plus_arg); hcross[k + (N/2+1)*(j + N*i)][s] *= cos(cross_arg); }
                //else if(s==1) { hplus[k + (N/2+1)*(j + N*i)][s] *= sin(plus_arg); hcross[k + (N/2+1)*(j + N*i)][s] *= sin(cross_arg); }

                //hplus[k + (N/2+1)*(j + N*i)][s] *= sqrt(2. * 4. * pow(m_bkgd_params.E, 2.));
                //hcross[k + (N/2+1)*(j + N*i)][s] *= sqrt(2. * 4. * pow(m_bkgd_params.E, 2.));
            }

            calc_transferse_vectors(i, j, k, N, mhat, nhat);
            for (int l=0; l<3; l++) for (int p=l; p<3; p++) for(int s=0; s<2; s++)
            {
                hk[lut[l][p]][k + (N/2+1)*(j + N*i)][s] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + (N/2+1)*(j + N*i)][s]
                                                    + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + (N/2+1)*(j + N*i)][s]) / sqrt(2.0);
            }
        }

        // If at a DC or Nyq point, enforce reality condition
        if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2)) // reality condition
        {
            hplus[k + (N/2+1)*(j + N*i)][1] = 0.;
            hcross[k + (N/2+1)*(j + N*i)][1] = 0.;
            for (int l=0; l<3; l++) for (int p=l; p<3; p++) { hk[lut[l][p]][k + (N/2+1)*(j + N*i)][1] = 0.; }
        }
    }

    //std::ofstream hkprint(m_params.print_path+"/h-k-printed.dat");
    //hkprint << std::fixed << setprecision(12);

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
    {
        apply_symmetry_rules(i, j, k, hplus, N);
        apply_symmetry_rules(i, j, k, hcross, N);
        for(int s=0; s<6; s++) { apply_symmetry_rules(i, j, k, hk[s], N); }

       /* for(int l=0; l<3; l++) for(int p=l; p<3; p++) for(int s=0; s<2; s++)
        {
            hkprint << hk[lut[l][p]][k + (N/2+1)*(j + N*i)][s] << ",";
            //hkprint << hplus[k + (N/2+1)*(j + N*i)][l] << "," << hcross[k + (N/2+1)*(j + N*i)][l] << ",";
        }
        hkprint << "\n";*/
    }

    //hkprint.close();
    //MayDay::Error("Printed file for comparison with stand-alone IC generator.");

    pout() << "Moving to configuration space.\n";

    fftw_execute(plan1);
    fftw_execute(plan2);
    for(int l=0; l<6; l++)
    { 
        fftw_execute(hij_plan[l]);
    }

    //std::ofstream hijprint(m_params.print_path+"/hij-printed.dat");
    //hijprint << std::fixed << setprecision(15);

    int Nc = m_params.N;
    int skip = (int)(N/Nc);
    double dx = m_params.L/Nc;

    std::vector<double> means(2, 0.);
    for(int i=0; i<Nc; i++) for(int j=0; j<Nc; j++) for(int k=0; k<Nc; k++)
    {
        /*if(m_spec_type == "position")
        {
            for(int l=0; l<3; l++) for(int p=l; p<3; p++)
            {
                hijprint << hx[lut[l][p]][k*skip + N * (j*skip + N * i*skip)] * norm << ",";
            }
            hijprint << "\n";
        }*/

        hplusx[(k + N * (j + N * i))*skip] *= norm;
        hcrossx[(k + N * (j + N * i))*skip] *= norm;

        means[0] += hplusx[(k + N * (j + N * i))*skip];
        means[1] += hcrossx[(k + N * (j + N * i))*skip];
    }
    //hijprint.close();
    //MayDay::Error("Check hij print file.");

    for(int s=0; s<2; s++) { means[s] /= pow(Nc, 3.); }

    std::vector<double> stdevs(2, 0.);
    for(int i=0; i<Nc; i++) for(int j=0; j<Nc; j++) for(int k=0; k<Nc; k++)
    {
        stdevs[0] += pow(hplusx[(k + N * (j + N * i))*skip] - means[0], 2.);
        stdevs[1] += pow(hcrossx[(k + N * (j + N * i))*skip] - means[1], 2.);
    }

    for(int s=0; s<2; s++)
    {
        stdevs[s] /= pow(Nc, 3.);
        stdevs[s] = sqrt(stdevs[s]);
    }

    if (m_spec_type == "position")
    {
    	ofstream pert_chars(m_params.print_path+"/IC-pert-level.dat");
	    if(!pert_chars) { MayDay::Error("Pert. IC characteristics file unopened."); }

    	pert_chars << "Mass scale [Mp]: " << m_bkgd_params.E << "\n";
    	pert_chars << "Length of box [Mp]: " << m_params.L * m_bkgd_params.E << "\n";
        pert_chars << "Full box resolution: " << N << "\n";
    	pert_chars << "Coarse box resolution: " << Nc << "\n";
    	pert_chars << "Amplitude: " << m_params.A << "\n";
        pert_chars << "Means: " << means[0] << ", " << means[1] << "\n";
    	pert_chars << "Std. deviation of plus pol. field: " << stdevs[0] << "\n";
    	pert_chars << "Std. deviation of cross pol. field: " << stdevs[1] << "\n";
    	pert_chars.close();
    }	
    stdevs.clear();


    // Free everything
    // (!!) note for c2r transforms, the original k-space array
    // is automatically destroyed by fftw_execute (annoying, I know...)
    free(hplusx);
    free(hcrossx);

    fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    for(int s=0; s<6; s++) { fftw_destroy_plan(hij_plan[s]); }

    pout() << "All memory but hx freed, starting BoxLoop\n";
    //if(m_spec_type == "velocity") { MayDay::Error("Check STDEVS print file."); }
}

inline int RandomField::flip_index(int I, int N) { return (int)abs(N - I); }
inline int RandomField::invert_index(int I, int N) { return (int)(N/2 - abs(N/2 - I)); }
inline int RandomField::invert_index_with_sign(int I, int N)
{
    if(I <= N/2) { return I; }
    else { return abs(N/2 - I) - N/2; }
}

inline void RandomField::apply_symmetry_rules(int i, int j, int k, double field[][2], int N)
{
    if (k==0 || k==N/2) 
    {
        if((i>N/2 && j==N/2) || (i==0 && j>N/2) || (i>N/2 && j==0) || (i==N/2 && j>N/2)) // Special lines
        {
            field[k + (N/2+1)*(j + N*i)][0] = field[k + (N/2+1)*(invert_index(j, N) + N*invert_index(i, N))][0];
            field[k + (N/2+1)*(j + N*i)][1] = -field[k + (N/2+1)*(invert_index(j, N) + N*invert_index(i, N))][1];
        }
        else if(j > N/2) // Special plane bulk
        {
            field[k + (N/2+1)*(j + N*i)][0] = field[k + (N/2+1)*(invert_index(j, N) + N*flip_index(i, N))][0];
            field[k + (N/2+1)*(j + N*i)][1] = -field[k + (N/2+1)*(invert_index(j, N) + N*flip_index(i, N))][1];
        }
    }
}

inline double RandomField::find_rayleigh_factor(double km, std::string spec_type)
{
    if(km < 1.e-12) { return 0.; } // P(k=0), for m=0

    double windowed_value = 0.;
    if (spec_type == "position")
    {
        windowed_value = (1.0/km/2.0 + (H0*H0/km/km/km)/2.0);
    }
    else if (m_spec_type == "velocity")
    {
        windowed_value = (km/2.0 - (H0*H0)/km/2.0 + H0*H0*H0*H0/km/km/km/2.0); 
    }

    // Apply the normalisation required to translate the scalar PS into tensor PS
    //windowed_value *= 2. * 4. * pow(m_bkgd_params.E, 2.); // 8/Mp where Mp is in units of the energy scale

    // Apply the tanh window function and the uniform draw
    windowed_value *= 0.5 * (1.0 - tanh(epsilon * (km - kstar)));
    return windowed_value;
}

inline void RandomField::calc_transferse_vectors(int x, int y, int z, int N, double MHat[3], double NHat[3], double a)
{
    double mh[3];
    double nh[3];
    for (int l=0; l<3; l++) { MHat[l] = 0.; NHat[l] = 0.; mh[l]=0.; nh[l]=0.;}

    if (a < 0 || a >= 2.*M_PI) 
    { 
        MayDay::Error("RandomField: Please choose a shift factor between 0 and 2 pi."); 
    }

    double X = x;
    double Y = y;
    double Z = z;

    if(x > N/2) { X = invert_index_with_sign(x, N); }
    if(y > N/2) { Y = invert_index_with_sign(y, N); }

    if (Z > 0.) 
    {
        if (X == 0. && Y == 0.) { mh[0] = 1.; mh[1] = 0.; mh[2] = 0.; 
                                  nh[0] = 0.; nh[1] = 1.; nh[2] = 0.; 
                                }

        else { mh[0] = Y/sqrt(X*X+Y*Y); mh[1] = -X/sqrt(X*X+Y*Y); mh[2] = 0.L;
               nh[0] = Z*X/sqrt(Z*Z*(X*X + Y*Y) + pow(X*X + Y*Y, 2.));
               nh[1] = Z*Y/sqrt(Z*Z*(X*X + Y*Y) + pow(X*X + Y*Y, 2.));
               nh[2] = -(X*X + Y*Y)/sqrt(Z*Z*(X*X + Y*Y) + pow(X*X + Y*Y, 2.)); 
             }
    }

    else if (abs(Y) > 0) { mh[0] = 0.; mh[1] = 0.; mh[2] = -1.;
                      nh[0] = -Y/sqrt(Y*Y + X*X);
                      nh[1] = X/sqrt(Y*Y + X*X);
                      nh[2] = 0.; 
                    }

    else if (abs(X) > 0) { mh[0] = 0.; mh[1] = 1.; mh[2] = 0.;
                      nh[0] = 0.; nh[1] = 0.; nh[2] = 1.;
                    }

    else if (X==0 && Y==0 && Z==0) { ; }

    else 
    {
        MayDay::Error("RandomField: Part of Fourier space is not covered by polarisation tensor calculation.");
    }

    for(int l=0; l<3; l++) { MHat[l] = cos(a)*mh[l] + sin(a)*nh[l]; NHat[l] = -sin(a)*mh[l] + cos(a)*nh[l]; }

    if (X != 0 && Y != 0 && Z != 0) { Test_norm(MHat); Test_norm(NHat); Test_orth(MHat, NHat); }
}

inline void RandomField::Test_norm(double vec[]) 
{
    double norm = 0;
    for (int i=0; i<3; i++) { norm += vec[i]*vec[i]; }

    if (abs(1. - norm) > 1.e-12) 
    {
        cout << "A basis vector is not normalised! Norm: " << norm << "\n";
        MayDay::Error();
    }
}

inline void RandomField::Test_orth(double vec1[], double vec2[]) 
{
    double orth = 0.;
    for(int i=0; i<3; i++) { orth += vec1[i]*vec2[i]; }

    if (abs(orth) > 1.e-12) 
    {
        cout << "Two basis vectors are not orthogonal! Orth factor: " << orth << "\n";
        MayDay::Error();
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
