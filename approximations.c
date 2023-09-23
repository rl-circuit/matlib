/*
    This library gives you the ability to handle C double approximation errors.
    Created by: Rares Luchian (lr.professional@outlook.com) in Anno Domini MMXXIII.
    This is FOSS. Enjoy! :))
*/

#define EPSILON 1e-8

void approxZero(double* num)
{
    if(-EPSILON <= *num && *num <= EPSILON)
        *num = 0;
}

