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
}

template <class data_t>
void RandomField::compute(Cell<data_t> current_cell) const
{
    Coordinates<data_t> coords(current_cell, m_params.L/m_params.N, m_params.center);

    // Pull out the grid parametersß
    int N = m_params.N;
    double L = m_params.L;
    double dx = L/N;

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
    if(i < 0)
    {
        i = N + i;
    }
    else if(i >= N)
    {
        i = i - N;
    }

    if(j < 0)
    {
        j = N + j;
    }
    else if(j >= N)
    {
        j = j - N;
    }

    if(k < 0)
    {
        k = N + k;
    }
    else if(k >= N)
    {
        k = k - N;
    }

    // The flattened position (leading with z?)
    int r = k + N*(j + N*i);

    if (r < 0)
    {
        cout << r << endl;
        MayDay::Error("Cell index value below zero.");
    }
    else if(r > pow(m_params.N, 3.))
    {
        cout << r << endl;
        MayDay::Error("Cell index greater than resolution^3 at coarsest level.");
    }

    if(r==1) {cout << "In compute: " << hx[0][r][0] << "," << m_params.A*hx[0][r][0] << "\n"; }

    if(std::isnan(hx[0][r][0])) { MayDay::Error("Values are nan in the RandomField compute"); }

    if(m_spec_type == "position")
    {
        //store tensor metric variables, g_ij = delta_ij + 1/2 h_ij
        current_cell.store_vars(1., c_h11);
        current_cell.store_vars(0., c_h12);
        current_cell.store_vars(0., c_h13);
        current_cell.store_vars(1., c_h22);
        current_cell.store_vars(0., c_h23);
        current_cell.store_vars(1., c_h33);
    }

    else if(m_spec_type == "velocity")
    {
        current_cell.store_vars(0., c_A11);
        current_cell.store_vars(0., c_A12);
        current_cell.store_vars(0., c_A13);
        current_cell.store_vars(0., c_A22);
        current_cell.store_vars(0., c_A23);
        current_cell.store_vars(0., c_A33);
    }


    /*if(m_spec_type == "position")
    {
        //store tensor metric variables, g_ij = delta_ij + 1/2 h_ij
        current_cell.store_vars(1. + m_params.A * 0.5 * hx[0][r][0], c_h11);
        current_cell.store_vars(m_params.A * 0.5 * hx[1][r][0], c_h12);
        current_cell.store_vars(m_params.A * 0.5 * hx[2][r][0], c_h13);
        current_cell.store_vars(1. + m_params.A * 0.5 * hx[4][r][0], c_h22);
        current_cell.store_vars(m_params.A * 0.5 * hx[5][r][0], c_h23);
        current_cell.store_vars(1. + m_params.A * 0.5 * hx[8][r][0], c_h33);
    }

    else if(m_spec_type == "velocity")
    {
        current_cell.store_vars(-m_params.A * 0.5 * hx[0][r][0], c_A11);
        current_cell.store_vars(-m_params.A * 0.5 * hx[1][r][0], c_A12);
        current_cell.store_vars(-m_params.A * 0.5 * hx[2][r][0], c_A13);
        current_cell.store_vars(-m_params.A * 0.5 * hx[4][r][0], c_A22);
        current_cell.store_vars(-m_params.A * 0.5 * hx[5][r][0], c_A23);
        current_cell.store_vars(-m_params.A * 0.5 * hx[8][r][0], c_A33);
    }*/

    else { MayDay::Error("Spec type entered is not a viable option."); }

    fftw_free(**hx);
}

//template <class data_t>
void RandomField::calc_spectrum()
{
    // Set up parameters based on L and N
    double kstar = 16.*(2.*M_PI/m_params.L);
    double epsilon = 2./m_params.L;
    double H0 = -3.0*sqrt((8.0 * M_PI/3.0/m_params.m_pl/m_params.m_pl)*(0.5*m_params.velocity*m_params.velocity + 0.5*pow(m_params.m * m_params.amplitude, 2.0)));
    double norm = pow(m_params.N, 3.);
    int N = 10;

    //double kstar = 0.4*pow((pow(N, 2.0) + pow(N, 2.0) + pow(N/2, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

    // Polarisation basis vectors
    double mhat[3] = {0., 0., 0.};
    double nhat[3] = {0., 0., 0.}; 

    // Allocate memory for hij and mode functions
    hk = (fftw_complex**) malloc(sizeof(fftw_complex*) * 9);
    hx = (fftw_complex**) malloc(sizeof(fftw_complex*) * 9);
    fftw_plan hij_plan[9];

    for(int l=0; l<9; l++)
    {
        hk[l] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
        hx[l] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
        hij_plan[l] = fftw_plan_dft_3d(N, N, N, hk[l], hx[l], FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    fftw_complex (*hplus);
    fftw_complex (*hcross);
    hplus = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    hcross = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);

    // Extra memory for reality check on h+
    fftw_complex (*hplusx);
    hplusx = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    fftw_plan plan1 = fftw_plan_dft_3d(N, N, N, hcross, hplusx, FFTW_BACKWARD, FFTW_ESTIMATE);
    

    // Set all arrays to 0
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++) for(int s=0; s<2; s++)
    {
        hplus[k + N*(j + N*i)][s] = 0.;
        hcross[k + N*(j + N*i)][s] = 0.;
        for (int m=0; m<9; m++)
        {
            hk[m][k + N*(j + N*i)][s] = 0.;
            hx[m][k + N*(j + N*i)][s] = 0.;
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

    // Set up random number generators
    int seed;
    if(m_spec_type == "position") { seed = 3539263; }
    else if(m_spec_type == "velocity") { seed = 7586572; }
    else { MayDay::Error("Field type incorrectly specified -- please enter position or velocity."); }

    default_random_engine engine(seed);
    uniform_real_distribution<double> theta_dist(0, 2*M_PI);
    uniform_real_distribution<double> sigma_dist(0, 1);

    double rayleigh_factors[2] = {0., 0.};
    double theta_factors[2] = {0., 0.};

    //cout << "Starting Loop 1 (k <= N/2).\n";

    for(int i=0; i<N/2+1; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        int J,K;
        if(k > N/2) { K = invert_index(k, N); }
        if(j > N/2) { J = invert_index(j, N); }

        double kmag = pow((pow(i, 2.0) + pow(J, 2.0) + pow(K, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        for(int s=0; s<2; s++) 
        { 
            theta_factors[s] = theta_dist(engine); 
            rayleigh_factors[s] = find_rayleigh_factor(kmag, kstar, epsilon, m_spec_type, H0, sigma_dist(engine));
        }

        //Start of with random numbers in all places
        hplus[k + N*(j + N*i)][0] = rayleigh_factors[0] * cos(theta_factors[0]);
        hplus[k + N*(j + N*i)][1] = rayleigh_factors[0] * sin(theta_factors[0]);

        hcross[k + N*(j + N*i)][0] = rayleigh_factors[1] * cos(theta_factors[1]);
        hcross[k + N*(j + N*i)][1] = rayleigh_factors[1] * sin(theta_factors[1]);

        calc_transferse_vectors(i, j, k, mhat, nhat);
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            hk[lut[l][p]][k + N*(j + N*i)][0] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][0]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][0]) / sqrt(2.0);

            hk[lut[l][p]][k + N*(j + N*i)][1] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][1]
                                                + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][1]) / sqrt(2.0);
        }
    }


    for(int i=0; i<N/2+1; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        // If you're at a special point 
        if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2)) 
        { 
            hplus[k + N*(j + N*i)][0] = rayleigh_factors[0] * cos(theta_factors[0]);
            hplus[k + N*(j + N*i)][1] = 0.;

            hcross[k + N*(j + N*i)][0] = rayleigh_factors[1] * cos(theta_factors[1]);
            hcross[k + N*(j + N*i)][1] = 0.;

            // Find hij from the mode functions
            calc_transferse_vectors(i, j, k, mhat, nhat);
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hk[lut[l][p]][k + N*(j + N*i)][0] = ((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][0]
                                                    + (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][0]) / sqrt(2.0);

                hk[lut[l][p]][k + N*(j + N*i)][1] = 0.;//((mhat[l]*mhat[p] - nhat[l]*nhat[p]) * hplus[k + N*(j + N*i)][1]
                                                    //+ (mhat[l]*nhat[p] + nhat[l]*mhat[p]) * hcross[k + N*(j + N*i)][1]) / sqrt(2.0);
            }
        }


        // If you're on a special axis or plane, excluing k > N/2
        else if(i==0 || i == N/2)
        {
            //cout << i << "," << j << "," << k << "\n";
            // Special axis on either k-norm plane
            if((j > N/2 && k == 0) || (j > N/2 && k == N/2) || (k > N/2 && j == 0) || (k > N/2 && j == N/2))
            {
                hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(invert_index(j, N) + N*i)][0];
                hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(invert_index(j, N) + N*i)][1];

                hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(invert_index(j, N) + N*i)][0];
                hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(invert_index(j, N) + N*i)][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(invert_index(j, N) + N*i)][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(invert_index(j, N) + N*i)][1];
                }

                //cout << i << "," << j << "," << k << "," << hplus[k + N*(j + N*i)][0] << "," << hplus[k + N*(j + N*i)][1] << "\n";
            }

            // Special plane on either k-norm plane
            else if(k > N/2)
            {
                hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
                hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];

                hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
                hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];
                }
            }
        }

        else if(j==0 || j == N/2)
        {
            // Special axis on either k-norm plane
            if((k > N/2 && i == 0) || (k > N/2 && i == N/2)) // (i > N/2 && k == 0) || (i > N/2 && k == N/2) ||
            {
                hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*invert_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*invert_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*invert_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*invert_index(i, N))][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(j + N*invert_index(i, N))][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(j + N*invert_index(i, N))][1];
                }
            }

            else if(k > N/2)
            {
                hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(j + N*flip_index(i, N))][1];
                }
            }
        }

        else if(k == 0 || k == N/2)
        {
            // Special axis on either k-norm plane
            if((i == 0 && j > N/2) || (i == N/2 && j > N/2)) // (i > N/2 && j == 0) || (i > N/2 || j == N/2)
            {
                hplus[k + N*(j + N*i)][0] = hplus[k + N*(invert_index(j, N) + N*invert_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[k + N*(invert_index(j, N) + N*invert_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[k + N*(invert_index(j, N) + N*invert_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[k + N*(invert_index(j, N) + N*invert_index(i, N))][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][k + N*(invert_index(j, N) + N*invert_index(i, N))][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][k + N*(invert_index(j, N) + N*invert_index(i, N))][1];
                }
            }
            // Special plane on either k-norm plane
            else if(j > N/2)
            {
                hplus[k + N*(j + N*i)][0] = hplus[k + N*(invert_index(j, N) + N*flip_index(i, N))][0];
                hplus[k + N*(j + N*i)][1] = -hplus[k + N*(invert_index(j, N) + N*flip_index(i, N))][1];

                hcross[k + N*(j + N*i)][0] = hcross[k + N*(invert_index(j, N) + N*flip_index(i, N))][0];
                hcross[k + N*(j + N*i)][1] = -hcross[k + N*(invert_index(j, N) + N*flip_index(i, N))][1];

                for(int s=0; s<9; s++)
                {
                    hk[s][k + N*(j + N*i)][0] = hk[s][k + N*(invert_index(j, N) + N*flip_index(i, N))][0];
                    hk[s][k + N*(j + N*i)][1] = -hk[s][k + N*(invert_index(j, N) + N*flip_index(i, N))][1];
                }
            }
        }
    }

    std::ofstream hpxcheck("./h-plus-mode-check.dat");

    for(int i=0; i<N/2+1; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        hpxcheck << i << "," << j << "," << k << ": ";
        hpxcheck << hplus[k + N*(j + N*i)][0] << "," << hplus[k + N*(j + N*i)][1] << "\n";
    }

    hpxcheck.close();

    MayDay::Error("Check file done");

    //cout << "Starting Loop 2 (k > N/2).\n";

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=N/2 + 1; k<N; k++)
    {
        int I,J;
        if(i > N/2) { I = flip_index(i, N); }
        if(j > N/2) { J = flip_index(j, N); }

        double kmag = pow((pow(I, 2.0) + pow(J, 2.0) + pow(k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        for(int s=0; s<2; s++) 
        { 
            theta_factors[s] = theta_dist(engine); 
            rayleigh_factors[s] = find_rayleigh_factor(kmag, kstar, epsilon, m_spec_type, H0, sigma_dist(engine));
        }

        if((j == 0 || j == N/2) && (i == 0 || i == N/2))
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*i)][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*i)][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*i)][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*i)][1];

            for(int s=0; s<9; s++)
            {
                hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(j + N*i)][0];
                hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(j + N*i)][1];
            }
        }

        /*else if(j == 0 || j == N/2)
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(j + N*flip_index(i, N))][1];

            for(int s=0; s<9; s++)
            {
                hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(j + N*flip_index(i, N))][0];
                hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(j + N*flip_index(i, N))][1];
            }
        }

        else if(i == 0 || i == N/2)
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];

            for(int s=0; s<9; s++)
            {
                hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*i)][0];
                hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*i)][1];
            }
        }*/

        else
        {
            hplus[k + N*(j + N*i)][0] = hplus[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][0];
            hplus[k + N*(j + N*i)][1] = -hplus[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][1];

            hcross[k + N*(j + N*i)][0] = hcross[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][0];
            hcross[k + N*(j + N*i)][1] = -hcross[invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][1];

            for(int s=0; s<9; s++)
            {
                hk[s][k + N*(j + N*i)][0] = hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][0];
                hk[s][k + N*(j + N*i)][1] = -hk[s][invert_index(k, N) + N*(flip_index(j, N) + N*flip_index(i, N))][1];
            }
        }
    }

    //cout << "Checking h+(x) for reality.\n";

    /*std::ofstream hpxcheck("./h-plus-im-check.dat");

    fftw_execute(plan1);
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        if(hplusx[k + N*(j + N*i)][1] > 1.e-12)
        {
            MayDay::Error("hx(x) is not yet real"); 
        }
        hpxcheck << i << "," << j << "," << k << ": " << hplusx[k + N*(j + N*i)][0] << "\n";
    }

    hpxcheck.close();*/

    for(int s=0; s<9; s++)
    {
        fftw_execute(hij_plan[s]);
    }

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++) for(int l=0; l<9; l++)
    {
        if(hx[l][k + N*(j + N*i)][1] > 1.e-12) 
        { 
            cout << i << "," << j << "," << k << ": " << hx[l][k + N*(j + N*i)][1] << "\n";
            MayDay::Error("hij(x) is not yet real"); 
        }
    }

    std::string name;
    if(m_spec_type == "position")
    {
        name = "./hk-position.dat";
    }
    else if(m_spec_type == "velocity")
    {
        name = "./hk-velocity.dat";
    }

    ofstream hx_check(name);

    int r;
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++)
    {
        r = k + N*(j + N*i);
        hx_check << hplus[k + N*(j + N*i)][0] << "," << hplus[k + N*(j + N*i)][1] << "," << hcross[k + N*(j + N*i)][0] << "," << hcross[k + N*(j + N*i)][1] << "\n";
    }

    //cout << "In calc: " << hx[0][1][0] << "\n";

    // Free everything
    fftw_free(*hplus);
    fftw_free(*hcross);
    //fftw_free(*hplusx);
    fftw_free(**hk);
    //fftw_free(**hx);

    fftw_destroy_plan(plan1);
    for(int s=0; s<9; s++)
    {
        fftw_destroy_plan(hij_plan[s]);
    }

    //cout << "Freed everything\n";
}

double RandomField::find_rayleigh_factor(double km, double ks, double ep, std::string spec_type, double H0, double uniform_draw)
{
    if(km - 0. < 1.e-12) { return 0.; }

    double windowed_value = 0.;
    if (spec_type == "position")
    {
        windowed_value = 0.5*km;//(0.5*(1.0/km + H0*H0/km/km/km));
    }
    else if (spec_type == "velocity")
    {
        windowed_value = (0.5*(km - H0*H0/km + H0*H0*H0*H0/km/km/km));
    }

    //windowed_value *= 0.5 * (1.0 - tanh(ep * (km - ks)));

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

    if (x != 0 && y != 0 && z != 0) 
    {
        Test_norm(MHat);
        Test_norm(NHat);
        Test_orth(MHat, NHat);
    }
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
