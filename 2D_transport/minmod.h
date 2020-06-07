
/* minmod limiter, eqn 16.52 p186 Numerical method for cons law leveques */
#ifndef MINMOD_H
inline __host__ __device__ double minmod(const double a, const double b)
{
    if ((abs(a) <= abs(b)) && (a * b >= 0)) {
        return a;
    } else if ((abs(a) > abs(b)) && (a * b >= 0)) {
        return b;
    } else {
        return 0;
    }
}
#define MINMOD_H
#endif
