#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "../src/doctest.h"
#include "../include/filter.h"
#include <numeric>

#ifdef DOCTEST_LIBRARY_INCLUDED



TEST_SUITE_BEGIN("filter:");


TEST_CASE("second order scalar filtering V4:"){

    using V = Vec4f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 1;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}



TEST_CASE("second order block filtering V4:"){

    using V = Vec4f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 131072;
    using T = float;
    constexpr static int M = 1;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}



TEST_CASE("4th order block filtering:"){

    using V = Vec8f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 2;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2,
                        xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}


TEST_CASE("second order multi-block filter V4:"){

    using V = Vec4f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 1;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}


TEST_CASE("second order iir filter V8:"){

    using V = Vec8f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 1;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}


TEST_CASE("second order iir filter V16:"){

    using V = Vec16f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 1;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}


TEST_CASE("4th order iir filter:"){

    using V = Vec8f;
    constexpr static int L = V::size();
    constexpr static int N = L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 2;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }

    T coefs[M][5] = {1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2,
                        xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}

TEST_CASE("8th order iir filter:"){

    using V = Vec8f;
    constexpr static int L = V::size();
    constexpr static int N = 4*L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 4;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }


    T coefs[M][5] = {1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}

TEST_CASE("16th order iir filter:"){

    using V = Vec8f;
    constexpr static int L = V::size();
    constexpr static int N = 4*L;
    constexpr static int vector_size = 16384;
    using T = float;
    constexpr static int M = 8;

    float b2 = 2, b1 = 1, a1 = 0.4, a2 = 0.5, xi2 = 1, xi1 = 2, yi1 = -3, yi2 = -5;

    std::array<T,vector_size> data{0},in,out;
    data[0] = 1; // pass an impulse response 

    std::array<T,vector_size> tmp,data_out = data;
    
    for (int round=0;round<M;round++){
        for (int n=0;n<vector_size;n++){

            if (n == 0) tmp[0] = data_out[0] + b2*xi2 + b1*xi1 + a2*yi2 + a1*yi1;
            else if (n == 1) tmp[1] = data_out[1] + b2*xi1 + b1*data_out[0] + a2*yi1 + a1*tmp[0];
            else tmp[n] = data_out[n] + b2*data_out[n-2] + b1*data_out[n-1] + a2*tmp[n-2] + a1*tmp[n-1];
        }
        data_out = tmp;
    }


    T coefs[M][5] = {1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2,
                     1,b1,b2,a1,a2}; 
    T inits[M][4] = {xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2,
                     xi1,xi2,yi1,yi2};

    auto _F = Filter<V,M,N>(coefs, inits);  

    in = data;

    _F(in.begin(),in.end(),out.begin());  
    
    for (auto r=0; r<vector_size; r++) CHECK(data_out[r] == doctest::Approx(out[r]));
}


TEST_SUITE_END();

#endif // doctest
