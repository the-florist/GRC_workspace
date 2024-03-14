/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(RANDOMFIELD_HPP_)
 #error "This file should only be included via RandomField.hpp"
 #endif

 #ifndef RANDOMFIELD_IMPL_HPP_
 #define RANDOMFIELD_IMPL_HPP_

 inline RandomField::RandomField(InitialScalarData::params_t a_params, std::string a_spec_type)
    : m_params(a_params), m_spec_type(a_spec_type)
{
    kstar = 8.*(2.*M_PI/m_params.L);
    epsilon = 0.05;
    H0 = 0.204692;//-3.0*sqrt((8.0 * M_PI/3.0/m_params.m_pl/m_params.m_pl)
            //*(0.5*m_params.velocity*m_params.velocity + 0.5*pow(m_params.m * m_params.amplitude, 2.0)));
    norm = pow(m_params.N, 3.);

    calc_spectrum();
}

template <class data_t>
void RandomField::compute(Cell<data_t> current_cell) const
{
    // Pull out the grid parameters
    int N = m_params.N;
    double L = m_params.L;
    double dx = L/N;

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
    if(i < 0) { i = N + i; }
    else if(i >= N) { i = i - N; }

    if(j < 0) { j = N + j; }
    else if(j >= N) { j = j - N; }

    if(k < 0) { k = N + k; }
    else if(k >= N) { k = k - N; }

    // The flattened position (leading with z)
    int r = k + N*(j + N*i);

    // Error trap, to make sure we stay in the domain of dependence
    if (r < 0)
    {
        cout << r << endl;
        MayDay::Error("RandomField: Cell index value below zero.");
    }
    else if(r > pow(m_params.N, 3.))
    {
        cout << r << endl;
        MayDay::Error("RandomField: Cell index greater than resolution^3 at coarsest level.");
    }

    if(m_spec_type == "position")
    {
        //store tensor metric variables, g_ij = delta_ij + 1/2 h_ij
        current_cell.store_vars(1. + m_params.A * hx[0][r], c_h11);
        current_cell.store_vars(m_params.A * hx[1][r], c_h12);
        current_cell.store_vars(m_params.A * hx[2][r], c_h13);
        current_cell.store_vars(1. + m_params.A * hx[4][r], c_h22);
        current_cell.store_vars(m_params.A * hx[5][r], c_h23);
        current_cell.store_vars(1. + m_params.A * hx[8][r], c_h33);
    }

    else if(m_spec_type == "velocity")
    {
        current_cell.store_vars(-m_params.A * hx[0][r], c_A11);
        current_cell.store_vars(-m_params.A * hx[1][r], c_A12);
        current_cell.store_vars(-m_params.A * hx[2][r], c_A13);
        current_cell.store_vars(-m_params.A * hx[4][r], c_A22);
        current_cell.store_vars(-m_params.A * hx[5][r], c_A23);
        current_cell.store_vars(-m_params.A * hx[8][r], c_A33);
    }

    else { MayDay::Error("RandomField: Spec type entered is not a viable option."); }

    // freeing the class memory that stores the config-space fields
    free(hx);
}

//template <class data_t>
void RandomField::calc_spectrum()
{
    int N = m_params.N;

    // Polarisation basis vectors and k
    double mhat[3] = {0., 0., 0.};
    double nhat[3] = {0., 0., 0.}; 
    double kmag = 0.;

    // Allocate memory for hij, mode functions and plans
    fftw_complex** hk;
    hk = (fftw_complex**) malloc(sizeof(fftw_complex*) * 9);
    hx = (double**) malloc(sizeof(double*) * 9);
    fftw_plan hij_plan[9];

    for(int l=0; l<9; l++)
    {
        hk[l] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * (N/2+1));
        hx[l] = (double*) malloc(sizeof(double) * N * N * N);
        hij_plan[l] = fftw_plan_dft_c2r_3d(N, N, N, hk[l], hx[l], FFTW_ESTIMATE);
    }

    fftw_complex (*hplus);
    fftw_complex (*hcross);
    hplus = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * (N/2+1));
    hcross = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * (N/2+1));

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

        for (int m=0; m<9; m++)
        {
            if(k <= N/2) 
            { 
                for(int s=0; s<2; s++) { hk[m][k + (N/2+1)*(j + N*i)][s] = 0.; }
            }
            hx[m][k + N*(j + N*i)] = 0.;
        }
    }

    // Setting the lut that maps polarisation vectors to 
    // polarisation tensors.
    int lut[3][3];  
    lut[0][0] = 0;
    lut[0][1] = 1;
    lut[0][2] = 2;
    lut[1][0] = 3;
    lut[1][1] = 4;
    lut[1][2] = 5;
    lut[2][0] = 6;
    lut[2][1] = 7;
    lut[2][2] = 8;

    // Set up random number generators (one independent seed per random draw)
    int seed;
    if(m_spec_type == "position") { seed = 3539263; }
    else if(m_spec_type == "velocity") { seed = 7586572; }
    else { MayDay::Error("RandomField: Please choose either 'position' or 'velocity' field type."); }

    default_random_engine engine(seed);
    uniform_real_distribution<double> theta_dist(0, 2*M_PI);
    uniform_real_distribution<double> sigma_dist(0, 1);

    cout << "Starting RandomField loop for " << m_spec_type << " field.\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
    {
        // Putting the 0 mode in the right spot in memory
        int I = i;
        int J = j;
        if(i > N/2) { I = invert_index(i, N); }
        if(j > N/2) { J = invert_index(j, N); }

        kmag = pow((pow(I, 2.0) + pow(J, 2.0) + pow(k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        // Start of with random numbers filling the entire array
        // Real parts of h+, hx and hij
        hplus[k + (N/2+1)*(j + N*i)][0] = find_rayleigh_factor(kmag, m_spec_type, sigma_dist(engine), 0) * cos(theta_dist(engine));
        hcross[k + (N/2+1)*(j + N*i)][0] = find_rayleigh_factor(kmag, m_spec_type, sigma_dist(engine), 0) * cos(theta_dist(engine));

        calc_transferse_vectors(i, j, k, mhat, nhat);
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            hk[lut[l][p]][k + (N/2+1)*(j + N*i)][0] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + (N/2+1)*(j + N*i)][0]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + (N/2+1)*(j + N*i)][0]) / sqrt(2.0);
        }

        // If at a DC or Nyq point, enforce reality condition
        if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2)) // reality condition
        {
            hplus[k + (N/2+1)*(j + N*i)][1] = 0.;
            hcross[k + (N/2+1)*(j + N*i)][1] = 0.;
            for (int l=0; l<3; l++) for (int p=0; p<3; p++) { hk[lut[l][p]][k + (N/2+1)*(j + N*i)][1] = 0.; }
        }
        // Else, fill the imaginary part of each field appropriately
        else
        {
            hplus[k + (N/2+1)*(j + N*i)][1] = find_rayleigh_factor(kmag, m_spec_type, sigma_dist(engine), 1) * sin(theta_dist(engine));
            hcross[k + (N/2+1)*(j + N*i)][1] = find_rayleigh_factor(kmag, m_spec_type, sigma_dist(engine), 1) * sin(theta_dist(engine));
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hk[lut[l][p]][k + (N/2+1)*(j + N*i)][1] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + (N/2+1)*(j + N*i)][1]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + (N/2+1)*(j + N*i)][1]) / sqrt(2.0);
            }
        }
    }

    cout << "All independent values have been assigned.\n Applying symmetry rules.\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
    {
        apply_symmetry_rules(i, j, k, hplus, N);
        apply_symmetry_rules(i, j, k, hcross, N);
        for(int s=0; s<9; s++) { apply_symmetry_rules(i, j, k, hk[s], N); }
    }

    cout << "Moving to configuration space.\n";

    fftw_execute(plan1);
    fftw_execute(plan2);
    for(int s=0; s<9; s++) { fftw_execute(hij_plan[s]); }

    // Free everything
    // (!!) note for c2r transforms, the original k-space array
    // is automatically destroyed by fftw_execute (annoying, I know...)
    free(hplusx);
    free(hcrossx);

    fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    for(int s=0; s<9; s++) { fftw_destroy_plan(hij_plan[s]); }
}

int RandomField::flip_index(int I, int N) { return (int)abs(N - I); }
int RandomField::invert_index(int I, int N) { return (int)(N/2 - abs(N/2 - I)); }

void RandomField::apply_symmetry_rules(int i, int j, int k, double field[][2], int N)
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

double RandomField::find_rayleigh_factor(double km, std::string spec_type, double uniform_draw, int comp)
{
    if(km < 1.e-12) { return 0.; } // P(k=0), for m=0

    double windowed_value = 0.;
    // See Mukanov-Sasaki mode function decomposition in: (forthcoming paper)
    if (spec_type == "position")
    {
        if(comp == 0) { windowed_value = (cos(km/H0) - H0 * sin(km/H0)/km)/sqrt(2. * km); }
        else if(comp == 1) { windowed_value = -(sin(km/H0) + H0 * cos(km/H0)/km)/sqrt(2. * km); }
        else { MayDay::Error("RandomField: component other than real or imaginary has been requested."); }
    }
    else if (spec_type == "velocity")
    {
        if(comp == 0) { windowed_value = (sin(km/H0) * (H0*H0/km - km) - H0 * cos(km/H0))/sqrt(2. * km); }
        else if(comp == 1) { windowed_value = (cos(km/H0) * (H0*H0/km - km) + H0 * sin(km/H0))/sqrt(2. * km); }
        else { MayDay::Error("RandomField: component other than real or imaginary has been requested."); }
    }

    // Apply the tanh window function and the uniform draw
    windowed_value *= 0.5 * (1.0 - tanh(epsilon * (km - kstar))) * sqrt(-2. * log(uniform_draw));
    return windowed_value;
}

void RandomField::calc_transferse_vectors(int x, int y, int z, double MHat[3], double NHat[3], double a)
{
    double mh[3];
    double nh[3];
    for (int l=0; l<3; l++) { MHat[l] = 0.; NHat[l] = 0.; mh[l]=0.; nh[l]=0.;}

    if (a < 0 || a >= 2.*M_PI) 
    { 
        MayDay::Error("RandomField: Please choose a shift factor between 0 and 2 pi."); 
    }

    if (z > 0.) 
    {
        if (x == 0. && y == 0.) { mh[0] = 1.; mh[1] = 0.; mh[2] = 0.; 
                                  nh[0] = 0.; nh[1] = 1.; nh[2] = 0.; 
                                }

        else { mh[0] = y/sqrt(x*x+y*y); mh[1] = -x/sqrt(x*x+y*y); mh[2] = 0.L;
               nh[0] = z*x/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.));
               nh[1] = z*y/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.));
               nh[2] = -(x*x + y*y)/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.)); 
             }
    }

    else if (y > 0) { mh[0] = 0.; mh[1] = 0.; mh[2] = -1.;
                      nh[0] = -y/sqrt(y*y + x*x);
                      nh[1] = x/sqrt(y*y + x*x);
                      nh[2] = 0.; 
                    }

    else if (x > 0) { mh[0] = 0.; mh[1] = 1.; mh[2] = 0.;
                      nh[0] = 0.; nh[1] = 0.; nh[2] = 1.;
                    }

    else if (x==0 && y==0 && z==0) { ; }

    else 
    {
        MayDay::Error("RandomField: Part of Fourier space is not covered by polarisation tensor calculation.");
    }

    for(int l=0; l<3; l++) { MHat[l] = cos(a)*mh[l] + sin(a)*nh[l]; NHat[l] = -sin(a)*mh[l] + cos(a)*nh[l]; }

    if (x != 0 && y != 0 && z != 0) { Test_norm(MHat); Test_norm(NHat); Test_orth(MHat, NHat); }
}

void RandomField::Test_norm(double vec[]) 
{
    double norm = 0;
    for (int i=0; i<3; i++) { norm += vec[i]*vec[i]; }

    if (abs(1. - norm) > 1.e-12) 
    {
        cout << "A basis vector is not normalised! Norm: " << norm << "\n";
        exit(EXIT_FAILURE);
    }
}

void RandomField::Test_orth(double vec1[], double vec2[]) 
{
    double orth = 0.;
    for(int i=0; i<3; i++) { orth += vec1[i]*vec2[i]; }

    if (abs(orth) > 1.e-12) 
    {
        cout << "Two basis vectors are not orthogonal! Orth factor: " << orth << "\n";
        exit(EXIT_FAILURE);
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
