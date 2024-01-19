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
    double kstar = 16.*(2.*M_PI/m_params.L);
    double epsilon = 2./m_params.L;
    double H0 = -3.0*sqrt((8.0 * M_PI/3.0/m_params.m_pl/m_params.m_pl)*(0.5*m_params.velocity*m_params.velocity + 0.5*pow(m_params.m * m_params.amplitude, 2.0)));
    double norm = pow(m_params.N, 3.);

    Coordinates<data_t> coords(current_cell, m_params.L/m_params.N, m_params.center);
}

/*double calc_h_real()
{
    return 0.;
}*/

void RandomField::calc_spectrum()
{
    int N = m_params.N;
    hk = (fftw_complex**) malloc(sizeof(fftw_complex*) * 9);
    for(int i=0; i<9; i++)
    {
            hk[i] = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    } 

    fftw_complex (*hplus);
    fftw_complex (*hcross);
    hplus = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);
    hcross = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N * N);

    int seed[2] = {7586572, 3539263};
    default_random_engine position_engine(seed[0]);
    default_random_engine velocity_engine(seed[1]);
    uniform_real_distribution<double> theta_engine(0, 2*M_PI);
    uniform_real_distribution<double> sigma_engine(0, 1);

    for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N/2+1; k++)
    {
        int I,J;
        if(i > N/2) { I = (int)(N/2 - abs(N/2 - i)); }
        if(j > N/2) { J = (int)(N/2 - abs(N/2 - j)); }

        double kmag = pow((pow(I, 2.0) + pow(J, 2.0) + pow(k, 2.0))*4.*M_PI*M_PI/m_params.L/m_params.L, 0.5);

        //if (i==0 || i == N/2) { hplus[k + N*(j + N*i)][0] = calc_h_real(); }
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/