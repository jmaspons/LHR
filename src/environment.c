#include <R.h>
#include <Rmath.h>
/*
 * SEASONAL EVIRONMENT follows a sinusoidal pattern.
 * @breedInterval in days
 */
int seasonalPattern(int events, double interval, double amplitude, double ** pattern)
{
    double * ptrPattern, day, mean;
    int i, interannual = 0;

    ptrPattern = Calloc(events, double);
    mean = 1 - amplitude / 2; // range(pattern) = [1 - amplitude, 1]

    for (i = 0, day = M_PI / 2; i < events; i++, day += interval * 2 * M_PI / 365) {
        ptrPattern[i] = sin(day) * amplitude + mean;

        if (day > 5 / 4 * M_PI)
            interannual = 1;
    }

    R_rsort(ptrPattern, events);
    * pattern = ptrPattern;

    return interannual? -events : events;
}

/*
seasonal.range<- function(min, max=1){
  mean<- (min + 1)/2
  amplitude<- max - min
// max-amplitude = min
// mean = max - amplitude/2
  return (data.frame(mean, amplitude))
}
*/
