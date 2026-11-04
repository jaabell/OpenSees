/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: process-global store for the load multiplier lambda produced
// by continuation methods (DisplacementControl, ArcLength, ...).
//
// WHY THIS EXISTS
// ---------------
// In a static analysis OpenSees historically stores lambda in
// Domain::currentTime: LoadControl does "t += dLambda", DisplacementControl
// does "t = lambda" where lambda solves the constraint equation. That makes
// the pseudo-time of a continuation stage a *solution variable*: it is not
// monotone (lambda decreases past a limit point, and may be negative), yet it
// is the value handed to every recorder by Domain::commit(). A model can
// therefore not be driven by amplitudes defined on a single global timeline.
//
// This store gives lambda its own slot, so Domain::currentTime can stay the
// one monotone pseudo-time. ContinuationTimeSeries reads this store; the
// continuation integrators write it.
//
// SEMANTICS
// ---------
// - a channel is "invalid" until a continuation integrator writes it. Reading
//   an invalid channel is a modelling error (a ContinuationTimeSeries used
//   without a continuation integrator) and yields a zero factor.
// - once written, the value PERSISTS after the stage ends. This is deliberate:
//   it is the propagation behaviour of a load applied by a previous stage,
//   the analogue of PathSeries' -useLast. It is not a stale-value bug.
//
// CAVEATS
// -------
// This is a per-process global, NOT part of the Domain state. It survives a
// database save/restore only because ContinuationTimeSeries carries the value
// of its channel through sendSelf/recvSelf and re-seeds the store on receive
// (never-written channels stay never-written). Under MPI every rank holds its
// own copy, so consistency relies on every rank's integrator computing the
// same lambda (which in turn relies on the distributed SOE returning the same
// solution vector on every rank).

#ifndef ContinuationLambda_h
#define ContinuationLambda_h

class OPS_ContinuationLambda
{
  public:
    enum { numChannels = 8 };

    // written by the continuation integrator, before applyLoadDomain()
    static void set(int channel, double value);

    // read by ContinuationTimeSeries::getFactor()
    static double get(int channel);

    // false until some continuation integrator has written this channel
    static bool isValid(int channel);

    // drop a channel back to the never-written state
    static void invalidate(int channel);
    static void invalidateAll(void);

    // identity of the integrator driving a channel, for diagnostics
    static void setOwner(int channel, int ownerTag);
    static int getOwner(int channel);

    static bool inRange(int channel);

  private:
    static double theValue[numChannels];
    static bool theValid[numChannels];
    static int theOwner[numChannels];
};

#endif
