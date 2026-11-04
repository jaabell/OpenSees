/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: see ContinuationTimeSeries.h

#include <ContinuationTimeSeries.h>
#include <ContinuationLambda.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>
#include <elementAPI.h>
#include <string.h>

void *
OPS_ContinuationTimeSeries(void)
{
  // timeSeries Continuation $tag <-channel $c> <-factor $f>
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING ContinuationTimeSeries - insufficient arguments\n";
    opserr << "      want: timeSeries Continuation tag? <-channel c?> <-factor f?>\n";
    return 0;
  }

  int tag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING ContinuationTimeSeries - failed to read tag\n";
    return 0;
  }

  int channel = 0;
  double cFactor = 1.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *arg = OPS_GetString();
    if (strcmp(arg, "-channel") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING ContinuationTimeSeries - -channel needs a value\n";
        return 0;
      }
      numData = 1;
      if (OPS_GetIntInput(&numData, &channel) < 0) {
        opserr << "WARNING ContinuationTimeSeries - failed to read channel\n";
        return 0;
      }
    } else if (strcmp(arg, "-factor") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING ContinuationTimeSeries - -factor needs a value\n";
        return 0;
      }
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &cFactor) < 0) {
        opserr << "WARNING ContinuationTimeSeries - failed to read factor\n";
        return 0;
      }
    }
    // unknown trailing args are ignored, as elsewhere in the series parsers
  }

  if (!OPS_ContinuationLambda::inRange(channel)) {
    opserr << "WARNING ContinuationTimeSeries - channel " << channel
           << " out of range [0," << (int)OPS_ContinuationLambda::numChannels - 1
           << "]\n";
    return 0;
  }

  return new ContinuationTimeSeries(tag, channel, cFactor);
}

ContinuationTimeSeries::ContinuationTimeSeries(int tag, int chan, double theFactor)
  : TimeSeries(tag, TSERIES_TAG_ContinuationTimeSeries),
    channel(chan), cFactor(theFactor), warned(false)
{

}

ContinuationTimeSeries::ContinuationTimeSeries()
  : TimeSeries(TSERIES_TAG_ContinuationTimeSeries),
    channel(0), cFactor(1.0), warned(false)
{

}

ContinuationTimeSeries::~ContinuationTimeSeries()
{

}

TimeSeries *
ContinuationTimeSeries::getCopy(void)
{
  return new ContinuationTimeSeries(this->getTag(), channel, cFactor);
}

double
ContinuationTimeSeries::getFactor(double pseudoTime)
{
  // pseudoTime is deliberately unused: the factor is the continuation
  // method's lambda, not a function of the analysis time.
  if (!OPS_ContinuationLambda::isValid(channel)) {
    if (warned == false) {
      opserr << "WARNING ContinuationTimeSeries (tag " << this->getTag()
             << ") - lambda channel " << channel << " has never been written "
             << "by a continuation integrator; returning a zero load factor.\n"
             << "      Use this series only on the reference load pattern of a "
             << "DisplacementControl/ArcLength stage.\n";
      warned = true;
    }
    return 0.0;
  }

  return cFactor * OPS_ContinuationLambda::get(channel);
}

int
ContinuationTimeSeries::sendSelf(int commitTag, Channel &theChannel)
{
  // The lambda value travels with the series so that a database save/restore
  // round-trips it: after a restore the process-global store is empty, and
  // this series is the only saved object that knows the channel. On the
  // parallel send path the value is redundant (the integrator's recvSelf
  // seeds the channel too) and recvSelf below never overwrites it.
  static Vector data(4);
  data(0) = (double)channel;
  data(1) = cFactor;
  data(2) = OPS_ContinuationLambda::isValid(channel) ? 1.0 : 0.0;
  data(3) = OPS_ContinuationLambda::get(channel); // zero if never written

  int dbTag = this->getDbTag();
  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "ContinuationTimeSeries::sendSelf() - channel failed to send data\n";
    return -1;
  }
  return 0;
}

int
ContinuationTimeSeries::recvSelf(int commitTag, Channel &theChannel,
                                 FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  int dbTag = this->getDbTag();
  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "ContinuationTimeSeries::recvSelf() - channel failed to receive data\n";
    channel = 0;
    cFactor = 1.0;
    return -1;
  }
  channel = (int)data(0);
  cFactor = data(1);
  warned = false;

  // seed the channel from the saved lambda, but never clobber a value a live
  // integrator has already written in this process (parallel send path).
  // When the script later re-creates the integrator, setContinuationTime
  // inherits this value, so a previous stage's load stays applied.
  if (data(2) == 1.0 && OPS_ContinuationLambda::isValid(channel) == false)
    OPS_ContinuationLambda::set(channel, data(3));

  return 0;
}

void
ContinuationTimeSeries::Print(OPS_Stream &s, int flag)
{
  s << "ContinuationTimeSeries - factor: " << cFactor
    << "  lambda channel: " << channel
    << "  current lambda: " << OPS_ContinuationLambda::get(channel)
    << (OPS_ContinuationLambda::isValid(channel) ? "" : " (channel never written)")
    << endln;
}
