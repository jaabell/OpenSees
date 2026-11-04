/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: a TimeSeries whose factor is the load multiplier lambda of a
// continuation method (DisplacementControl, ArcLength, ...) instead of a
// function of the pseudo-time.
//
// getFactor() DELIBERATELY IGNORES ITS ARGUMENT. It returns
//
//     cFactor * OPS_ContinuationLambda::get(channel)
//
// The continuation integrator writes that channel immediately before calling
// applyLoadDomain(), so the pattern carrying this series is scaled by lambda
// while Domain::currentTime is free to remain a monotone global timeline that
// every other pattern's amplitude is evaluated at.
//
// Put this series ONLY on the reference load pattern that the continuation
// method drives. Every other pattern (gravity, preloads, imposed
// displacements) keeps an ordinary amplitude on the global timeline.
//
// NAMING/DESIGN NOTE: a TimeSeries that does not depend on time is an abuse of
// the abstraction. It is used because TimeSeries::getFactor is the only
// polymorphic hook on the Domain -> LoadPattern -> load chain; the alternative
// is a new argument on a pure virtual shared by 17 classes, or a new member on
// LoadPattern. If that trade is ever regretted, LoadPattern owning lambda is
// the escape hatch.
//
// This series must not be used inside a GroundMotion: getDuration() and
// getTimeIncr() cannot be given a meaningful value.

#ifndef ContinuationTimeSeries_h
#define ContinuationTimeSeries_h

#include <TimeSeries.h>

class ContinuationTimeSeries : public TimeSeries
{
  public:
    ContinuationTimeSeries(int tag, int channel = 0, double cFactor = 1.0);
    ContinuationTimeSeries();
    ~ContinuationTimeSeries();

    TimeSeries *getCopy(void);

    // pseudoTime is ignored on purpose -- see the class comment
    double getFactor(double pseudoTime);

    double getDuration(void) { return 0.0; }
    double getPeakFactor(void) { return cFactor; }
    double getTimeIncr(double pseudoTime) { return 1.0; }

    int getChannel(void) const { return channel; }

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  private:
    int channel;      // which lambda channel this series reads
    double cFactor;   // scale applied on top of lambda
    bool warned;      // one-shot diagnostic for a never-written channel
};

#endif
