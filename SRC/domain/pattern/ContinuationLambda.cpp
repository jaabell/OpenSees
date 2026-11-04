/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: see ContinuationLambda.h

#include <ContinuationLambda.h>
#include <OPS_Globals.h>

double OPS_ContinuationLambda::theValue[OPS_ContinuationLambda::numChannels] = { 0.0 };
bool OPS_ContinuationLambda::theValid[OPS_ContinuationLambda::numChannels] = { false };
int OPS_ContinuationLambda::theOwner[OPS_ContinuationLambda::numChannels] = { -1 };

bool
OPS_ContinuationLambda::inRange(int channel)
{
  return (channel >= 0 && channel < numChannels);
}

void
OPS_ContinuationLambda::set(int channel, double value)
{
  if (!inRange(channel)) {
    opserr << "WARNING OPS_ContinuationLambda::set() - channel " << channel
           << " out of range [0," << (int)numChannels - 1 << "]\n";
    return;
  }
  theValue[channel] = value;
  theValid[channel] = true;
}

double
OPS_ContinuationLambda::get(int channel)
{
  if (!inRange(channel))
    return 0.0;
  return theValue[channel];
}

bool
OPS_ContinuationLambda::isValid(int channel)
{
  if (!inRange(channel))
    return false;
  return theValid[channel];
}

void
OPS_ContinuationLambda::invalidate(int channel)
{
  if (!inRange(channel))
    return;
  theValue[channel] = 0.0;
  theValid[channel] = false;
  theOwner[channel] = -1;
}

void
OPS_ContinuationLambda::invalidateAll(void)
{
  for (int i = 0; i < numChannels; i++)
    invalidate(i);
}

void
OPS_ContinuationLambda::setOwner(int channel, int ownerTag)
{
  if (!inRange(channel))
    return;
  theOwner[channel] = ownerTag;
}

int
OPS_ContinuationLambda::getOwner(int channel)
{
  if (!inRange(channel))
    return -1;
  return theOwner[channel];
}
