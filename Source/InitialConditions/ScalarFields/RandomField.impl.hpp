/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(RANDOMFIELD_HPP_)
 #error "This file should only be included via RandomField.hpp"
 #endif

 #ifndef RANDOMFIELD_IMPL_HPP_
 #define RANDOMFIELD_IMPL_HPP_

 inline RandomField::RandomField(InitialScalarData::params_t a_params)
    : m_params(a_params)
{
}

template <class data_t>
void RandomField::compute(Cell<data_t> current_cell) const
{
    Coordinates<data_t> coords(current_cell, m_params.L/m_params.N, m_params.center);
}

//template <class data_t>
void RandomField::calc_spectrum(std::string spec_type)
{
    // Set up parameters based on L and N
    double kstar = 16.*(2.*M_PI/m_params.L);
    double epsilon = 2./m_params.L;
    double H0 = -3.0*sqrt((8.0 * M_PI/3.0/m_params.m_pl/m_params.m_pl)*(0.5*m_params.velocity*m_params.velocity + 0.5*pow(m_params.m * m_params.amplitude, 2.0)));
    double norm = pow(m_params.N, 3.);
    int N = m_params.N;

    // Allocate memory for hij and mode functions
    hk = (fftw_complex**) malloc(sizeof(fftw_complex*) * 9);
    for(int i=0; i<9; i++)
    {
        hk[i] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    } 

    fftw_complex (*hplus);
    fftw_complex (*hcross);
    hplus = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    hcross = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);

    // Extra memory for reality check on h+
    fftw_complex (*hplusx);
    hplusx = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    fftw_plan plan1 = fftw_plan_dft_3d(N, N, N, hplus, hplusx, FFTW_BACKWARD, FFTW_ESTIMATE);

    /*fftw_plan (*hij_plan);
    hij_plan = (fftw_plan*) malloc(sizeof(fftw_plan) * 9);
    for(int l=0; l<9; l++)
    {
        hij_plan[l] = fftw_plan_dft_3d(N, N, N, hk[l], hx[l], FFTW_BACKWARD, FFTW_ESTIMATE);
    }*/

    // Set all arrays to 0
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++) for(int s=0; s<2; s++)
    {
        hplus[k + N*(j + N*i)][s] = 0.;
        hcross[k + N*(j + N*i)][s] = 0.;
        for (int m=0; m<9; m++)
        {
            hk[m][k + N*(j + N*i)][s] = 0.;
        }
    }

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

    // Set up random number generators
    int seed;
    if(spec_type == "position") { seed = 3539263; }
    else if(spec_type == "velocity") { seed = 7586572; }
    else { MayDay::Error("Field type incorrectly specified -- please enter position or velocity."); }

    default_random_engine engine(seed);
    uniform_real_distribution<double> theta_dist(0, 2*M_PI);
    uniform_real_distribution<double> sigma_dist(0, 1);

    double rayleigh_factors[2] = {0., 0.};
    double theta_factors[2] = {0., 0.};
    double mhat[3] = {0., 0., 0.};
    double nhat[3] = {0., 0., 0.}; 

    cout << "Starting Loop 1 (k <= N/2).\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N/2+1; k++)
    {
        int I,J;
        if(i > N/2) { I = invert_index(i, N); }
        if(j > N/2) { J = invert_index(j, N); }

        double kmag = pow((pow(I, 2.0) + pow(J, 2.0) + pow(k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        for(int s=0; s<2; s++) 
        { 
            theta_factors[s] = theta_dist(engine); 
            rayleigh_factors[s] = find_rayleigh_factor(kmag, kstar, epsilon, spec_type, H0, sigma_dist(engine));
        }

        // If you're at a special point 
        if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2)) 
        { 
            hplus[k + N*(j + N*i)][0] = rayleigh_factors[0] * cos(theta_factors[0]);
            hplus[k + N*(j + N*i)][1] = 0.;

            hcross[k + N*(j + N*i)][0] = rayleigh_factors[1] * cos(theta_factors[1]);
            hcross[k + N*(j + N*i)][1] = 0.;
        }

        // If you're on a special axis or plane, excluing k > N/2
        else if(k == 0 || k == N/2)
        {
            // Special axis on either k-norm plane
            if((i == 0 && j > N/2) || (i == N/2 && j > N/2) || (i > N/2 && j == 0) || (i > N/2 || j == N/2))
            {
                hplus[k + N*(j + N*i)][0] = hplus[k + N*(invert_index(j, N) + N*invert_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[k + N*(invert_index(j, N) + N*invert_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[k + N*(invert_index(j, N) + N*invert_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[k + N*(invert_index(j, N) + N*invert_index(i, N))][1];
            }
            // Special plane on either k-norm plane
            else if(j > N/2)
            {
                hplus[k + N*(j + N*i)][0] = hplus[k + N*(invert_index(j, N) + N*flip_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[k + N*(invert_index(j, N) + N*flip_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[k + N*(invert_index(j, N) + N*flip_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[k + N*(invert_index(j, N) + N*flip_index(i, N))][1];
            }
        }

        // The k < N/2 bulk
        else 
        {
            hplus[k + N*(j + N*i)][0] = rayleigh_factors[0] * cos(theta_factors[0]);
            hplus[k + N*(j + N*i)][1] = rayleigh_factors[0] * sin(theta_factors[0]);

            hcross[k + N*(j + N*i)][0] = rayleigh_factors[1] * cos(theta_factors[1]);
            hcross[k + N*(j + N*i)][1] = rayleigh_factors[1] * sin(theta_factors[1]);
        }

        // Find hij from the mode functions
        calc_transferse_vectors(i, j, k, mhat, nhat);
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            hk[lut[l][p]][k + N*(j + N*i)][0] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][0]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][0]) / sqrt(2.0);

            hk[lut[l][p]][k + N*(j + N*i)][1] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][1]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][1]) / sqrt(2.0);
        }
    }

    cout << "Starting Loop 2 (k > N/2).\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=N/2 + 1; k<N; k++)
    {
        int I,J;
        if(i > N/2) { I = flip_index(i, N); }
        if(j > N/2) { J = flip_index(j, N); }

        double kmag = pow((pow(I, 2.0) + pow(J, 2.0) + pow(k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        for(int s=0; s<2; s++) 
        { 
            theta_factors[s] = theta_dist(engine); 
            rayleigh_factors[s] = find_rayleigh_factor(kmag, kstar, epsilon, spec_type, H0, sigma_dist(engine));
        }

        if((j == 0 || j == N/2) && (i == 0 || i == N/2))
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*i)][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*i)][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*i)][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*i)][1];
        }
        else if(j == 0 || j == N/2)
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];
        }
        else if(i == 0 || i == N/2)
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];
        }
        else
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][1];
        }

        calc_transferse_vectors(i, j, k, mhat, nhat);
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            hk[lut[l][p]][k + N*(j + N*i)][0] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][0]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][0]) / sqrt(2.0);

            hk[lut[l][p]][k + N*(j + N*i)][1] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][1]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][1]) / sqrt(2.0);
        }
    }

    cout << "Checking h+(x) for reality.\n";

    fftw_execute(plan1);
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        if(hplusx[k + N*(j + N*i)][1] > 1e-12) { MayDay::Error("h+(x) is not yet real"); }
    }

    cout << "Freeing array\n";

    // Free everything
    fftw_free(*hplus);
    fftw_free(*hcross);
    //fftw_free(*hplusx);
    fftw_free(**hk);
    fftw_destroy_plan(plan1);
}

double RandomField::find_rayleigh_factor(double km, double ks, double ep, std::string spec_type, double H0, double uniform_draw)
{
    if(km - 0. < 1e-12) { return 0.; }

    double windowed_value = 0.;
    if (spec_type == "position")
    {
        windowed_value = (0.5*(1.0/km + H0*H0/km/km/km));
    }
    else if (spec_type == "velocity")
    {
        windowed_value = (0.5*(km - H0*H0/km + H0*H0*H0*H0/km/km/km));
    }

    windowed_value *= 0.5 * (1.0 - tanh(ep * (km - ks)));

    return windowed_value * sqrt(2./M_PI) * sqrt(-2. * log(uniform_draw));
}

int RandomField::flip_index(int I, int N) { return (int)abs(N - I); }
int RandomField::invert_index(int I, int N) { return (int)(N/2 - abs(N/2 - I)); }

void RandomField::calc_transferse_vectors(int x, int y, int z, double MHat[3], double NHat[3], double a)
{
    double mh[3];
    double nh[3];
    for (int l=0; l<3; l++) { MHat[l] = 0.; NHat[l] = 0.; mh[l]=0.; nh[l]=0.;}

    if (a < 0 || a >= 2.*M_PI) { MayDay::Error("Please choose a shift factor between 0 and 2 pi."); }

    if (z > 0.) 
    {
        if (x == 0. && y == 0.) 
        {
            mh[0] = 1.;
            mh[1] = 0.;
            mh[2] = 0.;

            nh[0] = 0.;
            nh[1] = 1.;
            nh[2] = 0.;
        }

        else 
        {
            mh[0] = y/sqrt(x*x+y*y);
            mh[1] = -x/sqrt(x*x+y*y);
            mh[2] = 0.L;

            nh[0] = z*x/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.));
            nh[1] = z*y/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.));
            nh[2] = -(x*x + y*y)/sqrt(z*z*(x*x + y*y) + pow(x*x + y*y, 2.));
        }

    }

    else if (y > 0) 
    {
        mh[0] = 0.;
        mh[1] = 0.;
        mh[2] = -1.;

        nh[0] = -y/sqrt(y*y + x*x);
        nh[1] = x/sqrt(y*y + x*x);
        nh[2] = 0.;

    }

    else if (x > 0) 
    {
        mh[0] = 0.;
        mh[1] = 1.;
        mh[2] = 0.;

        nh[0] = 0.;
        nh[1] = 0.;
        nh[2] = 1.;
    }

    else if (x==0 && y==0 && z==0) { ; }

    else 
    {
        MayDay::Error("Part of Fourier space is not covered by polarisation tensor calculation. Please examine your partition of this space.");
    }

    for(int l=0; l<3; l++) { MHat[l] = cos(a)*mh[l] + sin(a)*nh[l]; NHat[l] = -sin(a)*mh[l] + cos(a)*nh[l]; }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/