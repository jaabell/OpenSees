/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// Description: the interpreter commands of the contact subsystem, shared by
// the Python and Tcl interpreters (OpenSeesCommands, PythonWrapper,
// TclWrapper) and by the classic Tcl interpreter (SRC/tcl/commands.cpp):
//
//   contactSurface tag -slave nodeTag...
//   contactSurface tag -master nodesPerSeg nodeTag...
//   contactSurface tag -slave-segments nodesPerSeg nodeTag...
//   contact tag masterSurfTag slaveSurfTag <kn | kn kt mu | auto <kt mu>> <options>
//     options: -outward ox oy oz, -cell frac, -geomtan, -mu mu, -cohesion c,
//              -tauMax tau, -consistanttan, -visc muc
//     mortar:  -mortar, -epsN auto|val, -augTol tol, -maxAug n, -ngp order,
//              -tie, -epsTie auto|val, -epsT auto|val (friction: -mu, -cohesion,
//              -tauMax, -consistanttan as above)
//              -augment commit|request|never, -maxGap d, -gapOffset g0,
//              -adjust <tol>, -smoothN g0, -smoothT r
//     edge-to-edge (with -mortar): -edgeedge, -edgeKn auto|val, -edgeBand d,
//              -edgeMu mu, -edgeKt kt, -edgeCohesion c, -edgeTauMax tau,
//              -edgeConsistentTan, -edgeAlm, -edgeAugTol tol
//     -consistentTan is accepted as an alias of -consistanttan.
//   contactPlane tag slaveSurfTag nx ny nz px py pz kn <-visc muc>
//   ladrunoContactInfo            -> numContacts numCommits numReverts numMortarContacts
//   ladrunoContactForce slaveTag  -> summed normal contact force on the node
//   ladrunoMortarPenetration      -> max mortar penetration on the active set
//   ladrunoMortarTieResidual      -> max weighted relative displacement of the mortar ties
//   ladrunoEdgePenetration        -> max penetration of the -edgeAlm edge pairs in contact
//   ladrunoBeginAugment, ladrunoEndAugment
//                                 -> open / close a held-load augmentation sweep
//
// The definitions are stored on the LadrunoContactDomain of the current
// Domain, created on first use; `constraints LadrunoContact` builds the
// contact adapters from them.

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <ID.h>
#include <cstring>
#include "LadrunoContactDomain.h"
#include "LadrunoContactSurface.h"

// get the contact engine of the current Domain, creating it on first use
static LadrunoContactDomain *ladrunoGetOrCreateContactDomain()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        cd = new LadrunoContactDomain();
        theDomain->setLadrunoContactDomain(cd);
    }
    return cd;
}

// Count the numeric tokens at the current parse position, up to cap, without
// consuming them. A token counts as an option flag only when '-' is followed
// by a letter, so a negative number is a number. OPS_GetString() consumes the
// token on both interpreters and, on the Python interpreter, returns a fixed
// non-flag string for a numeric argument, so the loop is bounded by the number
// of remaining arguments and the cursor is rewound by the tokens read.
static int ladrunoCountLeadingNumbers(int cap)
{
    int nrem = OPS_GetNumRemainingInputArgs();
    if (cap > nrem) cap = nrem;
    int consumed = 0, nnum = 0;
    bool hitFlag = false;
    while (consumed < cap && !hitFlag) {
        const char *tok = OPS_GetString();
        consumed++;
        bool isFlag = (tok != 0 && tok[0] == '-' &&
                       ((tok[1] >= 'a' && tok[1] <= 'z') ||
                        (tok[1] >= 'A' && tok[1] <= 'Z')));
        if (isFlag) hitFlag = true;
        else        nnum++;
    }
    if (consumed > 0)
        OPS_ResetCurrentInputArg(-consumed);      // un-read every peeked token
    return nnum;
}

// contactSurface tag (-slave | -master nodesPerSeg | -slave-segments nodesPerSeg) nodeTag...
int OPS_LadrunoContactSurface()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - contactSurface tag (-slave | -master nodesPerSeg | "
                  "-slave-segments nodesPerSeg) nodeTag...\n";
        return -1;
    }
    int tag, n = 1;
    if (OPS_GetIntInput(&n, &tag) < 0) {
        opserr << "WARNING contactSurface - could not read tag\n";
        return -1;
    }
    const char *kindStr = OPS_GetString();
    LadrunoContactSurface::Kind kind;
    int nodesPerSeg = 0;
    if (strcmp(kindStr, "-slave") == 0) {
        kind = LadrunoContactSurface::SLAVE_NODES;
    } else if (strcmp(kindStr, "-master") == 0) {
        kind = LadrunoContactSurface::MASTER_SEGMENTS;
        n = 1;
        if (OPS_GetIntInput(&n, &nodesPerSeg) < 0 || nodesPerSeg < 3 || nodesPerSeg > 4) {
            opserr << "WARNING contactSurface -master - need nodesPerSeg (3 for a triangle, "
                      "4 for a quadrilateral facet)\n";
            return -1;
        }
    } else if (strcmp(kindStr, "-slave-segments") == 0) {
        // a faceted slave surface, for mortar contact
        kind = LadrunoContactSurface::SLAVE_SEGMENTS;
        n = 1;
        if (OPS_GetIntInput(&n, &nodesPerSeg) < 0 || nodesPerSeg < 3 || nodesPerSeg > 4) {
            opserr << "WARNING contactSurface -slave-segments - need nodesPerSeg (3 for a "
                      "triangle, 4 for a quadrilateral facet)\n";
            return -1;
        }
    } else {
        opserr << "WARNING contactSurface - kind must be -slave, -master or -slave-segments\n";
        return -1;
    }
    int nrem = OPS_GetNumRemainingInputArgs();
    if (nrem <= 0) {
        opserr << "WARNING contactSurface - no node tags given\n";
        return -1;
    }
    ID nodeTags(nrem);
    for (int i = 0; i < nrem; i++) {
        int v, one = 1;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING contactSurface - could not read node tag\n";
            return -1;
        }
        nodeTags(i) = v;
    }
    // Every node must exist when the surface is declared. A node missing at
    // analysis time can then only have been removed by the analysis, which
    // the handler prunes; a typo fails here, at the line that contains it.
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) {
        opserr << "WARNING contactSurface - no Domain\n";
        return -1;
    }
    for (int i = 0; i < nodeTags.Size(); i++) {
        Node *nd = theDomain->getNode(nodeTags(i));
        if (nd == 0) {
            opserr << "WARNING contactSurface " << tag << " - node " << nodeTags(i)
                   << " does not exist. Define every node before the surface.\n";
            return -1;
        }
        if (nd->getCrds().Size() != 3) {
            opserr << "WARNING contactSurface " << tag << " - node " << nodeTags(i)
                   << " has " << nd->getCrds().Size() << " coordinates; contact surfaces "
                      "need 3D nodes\n";
            return -1;
        }
    }

    LadrunoContactDomain *cd = ladrunoGetOrCreateContactDomain();
    if (cd == 0) return -1;
    LadrunoContactSurface *s = new LadrunoContactSurface(tag, kind, nodeTags, nodesPerSeg);
    if (cd->addSurface(s) < 0) { delete s; return -1; }
    return 0;
}

// contact tag masterSurfTag slaveSurfTag <kn | kn kt mu | auto <kt mu>>
//         <-outward ox oy oz> <-cell frac> <-geomtan>
//         <-mu mu> <-cohesion c> <-tauMax tau> <-consistanttan> <-visc muc>
int OPS_LadrunoContact()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - contact tag masterSurfTag slaveSurfTag <kn kt mu> <-outward ox oy oz>\n";
        return -1;
    }
    int idata[3], n = 3;
    if (OPS_GetIntInput(&n, idata) < 0) {
        opserr << "WARNING contact - could not read tag/master/slave\n";
        return -1;
    }
    double kn = 0.0, kt = 0.0, mu = 0.0;
    bool knAuto = false;
    // The kn slot takes a number (`kn` or `kn kt mu`) or the word `auto`
    // (optionally followed by `kt mu`). It may be omitted; the next token is
    // then an option flag or the end of the arguments.
    if (OPS_GetNumRemainingInputArgs() >= 1) {
        const char *peek = OPS_GetString();   // consume and classify
        bool isFlag = (peek != 0 && peek[0] == '-');
        bool isAuto = (peek != 0 && strcmp(peek, "auto") == 0);
        if (isAuto) {
            knAuto = true;                    // kn is sized at handle()
            if (OPS_GetNumRemainingInputArgs() >= 2) {
                const char *p2 = OPS_GetString();
                bool flag2 = (p2 != 0 && p2[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!flag2) {
                    double d[2]; int m = 2;
                    if (OPS_GetDoubleInput(&m, d) < 0) {
                        opserr << "WARNING contact - could not read kt mu after 'auto'\n";
                        return -1;
                    }
                    kt = d[0]; mu = d[1];
                }
            }
        } else if (isFlag) {
            OPS_ResetCurrentInputArg(-1);     // un-read it for the option loop
        } else {
            OPS_ResetCurrentInputArg(-1);     // un-read it and read the numbers
            int m = ladrunoCountLeadingNumbers(3);
            if (m == 2) {
                opserr << "WARNING contact - give kn alone or kn kt mu\n";
                return -1;
            }
            double d[3] = {0.0, 0.0, 0.0};
            if (OPS_GetDoubleInput(&m, d) < 0) {
                opserr << "WARNING contact - could not read kn kt mu\n";
                return -1;
            }
            kn = d[0]; kt = d[1]; mu = d[2];
        }
    }
    bool hasOutward = false;
    double outward[3] = {0.0, 0.0, 0.0};
    double cellFrac = 1.0;
    bool consistentNormal = false;
    bool consistentTan = false;
    double cohesion = 0.0, tauMax = 0.0;
    double muc = 0.0;
    bool hasCell = false;
    // mortar segment-to-segment contact
    bool isMortar = false;          // -mortar
    bool hasMortarOpt = false;      // a mortar-only option was given
    double epsN = 0.0;              // -epsN value
    bool epsNAuto = false;          // -epsN auto
    double augTol = 1e-8;           // -augTol
    int maxAug = 20, ngp = 2;       // -maxAug, -ngp
    bool isTie = false;             // -tie
    double epsT = 0.0;              // -epsT value
    bool epsTAuto = false;          // -epsT auto
    // controls of the mortar contact (setMortarContactOptions)
    int    augmentMode = LadrunoContactDomain::AUG_COMMIT;
    bool   hasControls = false;
    double maxGap = 0.0, gapOffset = 0.0, adjustTol = 0.0;
    bool   adjust = false;
    double smoothN = 0.0, smoothT = 0.0;   // smoothed law (setMortarSmoothing)
    // edge-to-edge contact (setMortarEdgeEdge)
    bool   edgeEdge = false, edgeKnAuto = false, edgeConsistentTan = false, edgeAlm = false;
    bool   hasEdgeOpt = false;     // an -edge* option other than -edgeedge was given
    double edgeKn = 0.0, edgeBand = 0.0, edgeMu = 0.0, edgeKt = 0.0;
    double edgeCohesion = 0.0, edgeTauMax = 0.0, edgeAugTol = 0.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (opt != 0 && strcmp(opt, "-mortar") == 0) {
            isMortar = true;
        } else if (opt != 0 && strcmp(opt, "-edgeedge") == 0) {
            // edge-to-edge contact for facets that meet edge-on (mortar only)
            edgeEdge = true;
        } else if (opt != 0 && strcmp(opt, "-edgeKn") == 0) {
            // edge penalty; auto sizes it from the solid that owns the master facet
            const char *e = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
            if (e != 0 && strcmp(e, "auto") == 0) {
                edgeKnAuto = true;
            } else {
                if (e != 0) OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (e == 0 || OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -edgeKn - need auto or a value\n";
                    return -1;
                }
                edgeKn = v[0];
            }
            hasEdgeOpt = true;
        } else if (opt != 0 && (strcmp(opt, "-edgeBand") == 0 || strcmp(opt, "-edgeMu") == 0 ||
                                strcmp(opt, "-edgeKt") == 0 || strcmp(opt, "-edgeCohesion") == 0 ||
                                strcmp(opt, "-edgeTauMax") == 0 ||
                                strcmp(opt, "-edgeAugTol") == 0)) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact " << opt << " - need a value\n";
                return -1;
            }
            if (strcmp(opt, "-edgeBand") == 0)          edgeBand = v[0];
            else if (strcmp(opt, "-edgeMu") == 0)       edgeMu = v[0];
            else if (strcmp(opt, "-edgeKt") == 0)       edgeKt = v[0];
            else if (strcmp(opt, "-edgeCohesion") == 0) edgeCohesion = v[0];
            else if (strcmp(opt, "-edgeTauMax") == 0)   edgeTauMax = v[0];
            else                                        edgeAugTol = v[0];
            hasEdgeOpt = true;
        } else if (opt != 0 && strcmp(opt, "-edgeConsistentTan") == 0) {
            edgeConsistentTan = true;
            hasEdgeOpt = true;
        } else if (opt != 0 && strcmp(opt, "-edgeAlm") == 0) {
            // augmented Lagrangian on the edge pairs (one multiplier per pair)
            edgeAlm = true;
            hasEdgeOpt = true;
        } else if (opt != 0 && strcmp(opt, "-edgeSoft") == 0) {
            opserr << "WARNING contact -edgeSoft - the soft penalty is sized from the time step "
                      "of an explicit integrator, and no integrator in this build reports one; "
                      "use -edgeKn\n";
            return -1;
        } else if (opt != 0 && strcmp(opt, "-tie") == 0) {
            // a permanent mesh tie (mortar only, checked below)
            isTie = true;
        } else if (opt != 0 && (strcmp(opt, "-epsN") == 0 || strcmp(opt, "-epsTie") == 0)) {
            // normal penalty of the mortar augmented Lagrangian, or the penalty
            // of a tie (-epsTie, the same slot: a tie has one penalty); auto
            // sizes it from the solid that owns the master facet, as kn auto
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact " << opt << " - need auto or a value\n";
                return -1;
            }
            const char *e = OPS_GetString();
            if (e != 0 && strcmp(e, "auto") == 0) {
                epsNAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact " << opt << " - need auto or a value\n";
                    return -1;
                }
                epsN = v[0];
            }
            hasMortarOpt = true;
        } else if (opt != 0 && strcmp(opt, "-augment") == 0) {
            // when the Uzawa updates of the mortar multipliers run
            const char *m = (OPS_GetNumRemainingInputArgs() > 0) ? OPS_GetString() : 0;
            if (m != 0 && strcmp(m, "commit") == 0)       augmentMode = LadrunoContactDomain::AUG_COMMIT;
            else if (m != 0 && strcmp(m, "request") == 0) augmentMode = LadrunoContactDomain::AUG_REQUEST;
            else if (m != 0 && strcmp(m, "never") == 0)   augmentMode = LadrunoContactDomain::AUG_NEVER;
            else {
                opserr << "WARNING contact -augment - need commit|request|never\n";
                return -1;
            }
            hasControls = true;
        } else if (opt != 0 && (strcmp(opt, "-smoothN") == 0 || strcmp(opt, "-smoothT") == 0)) {
            // band g0 > 0 of the smoothed normal onset / ratio 0 < r < 1 of the
            // rounded stick/slip corner
            bool isN = (strcmp(opt, "-smoothN") == 0);
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0 || v[0] <= 0.0 || (!isN && v[0] >= 1.0)) {
                opserr << (isN ? "WARNING contact -smoothN - need a band g0 > 0 (a length)\n"
                               : "WARNING contact -smoothT - need a ratio 0 < r < 1\n");
                return -1;
            }
            if (isN) smoothN = v[0];
            else     smoothT = v[0];
        } else if (opt != 0 && strcmp(opt, "-maxGap") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0 || v[0] <= 0.0) {
                opserr << "WARNING contact -maxGap - need a distance > 0\n";
                return -1;
            }
            maxGap = v[0];
            hasControls = true;
        } else if (opt != 0 && strcmp(opt, "-gapOffset") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -gapOffset - need a value (negative = interference)\n";
                return -1;
            }
            gapOffset = v[0];
            hasControls = true;
        } else if (opt != 0 && strcmp(opt, "-adjust") == 0) {
            // optional tolerance; a flag or the end of the arguments: adjust all
            adjust = true;
            hasControls = true;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *p = OPS_GetString();
                bool isFlag = (p != 0 && p[0] == '-' &&
                               ((p[1] >= 'a' && p[1] <= 'z') || (p[1] >= 'A' && p[1] <= 'Z')));
                OPS_ResetCurrentInputArg(-1);
                if (!isFlag) {
                    double v[1]; int m = 1;
                    if (OPS_GetDoubleInput(&m, v) < 0 || v[0] <= 0.0) {
                        opserr << "WARNING contact -adjust - optional tolerance must be a value "
                                  "> 0 (or omit it to adjust every paired node)\n";
                        return -1;
                    }
                    adjustTol = v[0];
                }
            }
        } else if (opt != 0 && strcmp(opt, "-epsT") == 0) {
            // tangential penalty of the mortar friction; auto: the normal penalty
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact -epsT - need auto or a value\n";
                return -1;
            }
            const char *e = OPS_GetString();
            if (e != 0 && strcmp(e, "auto") == 0) {
                epsTAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -epsT - need auto or a value\n";
                    return -1;
                }
                epsT = v[0];
            }
            hasMortarOpt = true;
        } else if (opt != 0 && strcmp(opt, "-augTol") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -augTol - need a value\n";
                return -1;
            }
            augTol = v[0];
            hasMortarOpt = true;
        } else if (opt != 0 && (strcmp(opt, "-maxAug") == 0 || strcmp(opt, "-ngp") == 0)) {
            int v[1]; int m = 1;
            if (OPS_GetIntInput(&m, v) < 0) {
                opserr << "WARNING contact " << opt << " - need an integer\n";
                return -1;
            }
            if (strcmp(opt, "-maxAug") == 0) maxAug = v[0];
            else                             ngp = v[0];
            hasMortarOpt = true;
        } else if (opt != 0 && strcmp(opt, "-outward") == 0) {
            // direction toward the slave's allowed half-space
            double o[3]; int m = 3;
            if (OPS_GetDoubleInput(&m, o) < 0) {
                opserr << "WARNING contact -outward - need ox oy oz\n";
                return -1;
            }
            outward[0] = o[0]; outward[1] = o[1]; outward[2] = o[2];
            hasOutward = true;
        } else if (opt != 0 && strcmp(opt, "-cell") == 0) {
            // broad-phase cell size as a multiple of the median segment diagonal;
            // a huge value gives one bucket, i.e. brute force
            double f[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, f) < 0) {
                opserr << "WARNING contact -cell - need a positive frac\n";
                return -1;
            }
            cellFrac = f[0];
            hasCell = true;
        } else if (opt != 0 && strcmp(opt, "-geomtan") == 0) {
            // geometric part of the consistent normal tangent, for curved or
            // large-sliding masters; symmetric, so any solver works
            consistentNormal = true;
        } else if (opt != 0 && (strcmp(opt, "-mu") == 0 || strcmp(opt, "-cohesion") == 0 ||
                                strcmp(opt, "-tauMax") == 0)) {
            // friction cone cap = min(mu N + cohesion, tauMax); -mu replaces the
            // positional mu
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact " << opt << " - need a value\n";
                return -1;
            }
            if (strcmp(opt, "-mu") == 0)            mu = v[0];
            else if (strcmp(opt, "-cohesion") == 0) cohesion = v[0];
            else                                    tauMax = v[0];
        } else if (opt != 0 && (strcmp(opt, "-consistanttan") == 0 ||
                                strcmp(opt, "-consistentTan") == 0)) {
            // the non-symmetric consistent friction tangent (quadratic Newton
            // convergence); a symmetric solver reads one triangle only and
            // silently drops the coupling, so a non-symmetric one is needed
            consistentTan = true;
            opserr << "WARNING contact -consistanttan: the non-symmetric consistent "
                      "friction tangent needs a non-symmetric solver (system FullGeneral/"
                      "UmfPack/BandGeneral); symmetric solvers will silently corrupt it.\n";
        } else if (opt != 0 && strcmp(opt, "-visc") == 0) {
            // viscous normal stabilization: a dashpot muc on the normal gap rate
            // while in contact, damping chatter in impact and rocking
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -visc - need a coefficient value\n";
                return -1;
            }
            muc = v[0];
        } else if (opt != 0 && strcmp(opt, "-soft") == 0) {
            opserr << "WARNING contact -soft - the soft penalty is sized from the time step "
                      "of an explicit integrator, and no integrator in this build reports "
                      "one; use kn or auto\n";
            return -1;
        } else {
            opserr << "WARNING contact - unexpected token '" << (opt ? opt : "")
                   << "' (expected -outward, -cell, -geomtan, -mu, -cohesion, -tauMax, "
                      "-consistanttan, -visc, -mortar, -tie, -epsN, -epsTie, -epsT, -augTol, "
                      "-maxAug, -ngp, -augment, -maxGap, -gapOffset, -adjust, -smoothN, "
                      "-smoothT, -edgeedge, -edgeKn, -edgeBand, -edgeMu, -edgeKt, "
                      "-edgeCohesion, -edgeTauMax, -edgeConsistentTan, -edgeAlm, -edgeAugTol "
                      "or end of arguments)\n";
            return -1;
        }
    }
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        opserr << "WARNING contact - define a contactSurface first\n";
        return -1;
    }
    if (hasControls && !isMortar) {
        opserr << "WARNING contact -augment/-maxGap/-gapOffset/-adjust are -mortar options; "
                  "add -mortar or remove them\n";
        return -1;
    }
    if ((smoothN > 0.0 || smoothT > 0.0) && !isMortar) {
        opserr << "WARNING contact -smoothN/-smoothT are -mortar options; add -mortar or "
                  "remove them\n";
        return -1;
    }
    if (isTie && !isMortar) {
        opserr << "WARNING contact -tie requires -mortar (mesh tying is a mortar formulation)\n";
        return -1;
    }
    if (edgeEdge && !isMortar) {
        opserr << "WARNING contact -edgeedge requires -mortar (edge-to-edge contact is routed "
                  "from the facet pairs of a mortar contact)\n";
        return -1;
    }
    if (hasEdgeOpt && !edgeEdge) {
        opserr << "WARNING contact -edgeKn/-edgeBand/-edgeMu/-edgeKt/-edgeCohesion/-edgeTauMax/"
                  "-edgeConsistentTan/-edgeAlm/-edgeAugTol need -edgeedge\n";
        return -1;
    }
    if (edgeEdge && isTie) {
        opserr << "WARNING contact -edgeedge does not apply to -tie (a tie bonds the facets it "
                  "overlaps)\n";
        return -1;
    }
    if (edgeEdge && (edgeKn < 0.0 || edgeBand < 0.0 || edgeMu < 0.0 || edgeKt < 0.0 ||
                     edgeCohesion < 0.0)) {
        opserr << "WARNING contact -edgeKn, -edgeBand, -edgeMu, -edgeKt and -edgeCohesion must "
                  "be >= 0\n";
        return -1;
    }
    if (edgeEdge && edgeTauMax > 0.0 && edgeMu <= 0.0 && edgeCohesion <= 0.0) {
        opserr << "WARNING contact -edgeTauMax without -edgeMu or -edgeCohesion gives the cone "
                  "min(mu*N + c, tauMax) = 0 (free slip); for a shear-capped bond use "
                  "-edgeCohesion <c>, optionally with -edgeTauMax\n";
        return -1;
    }
    if (edgeEdge && edgeConsistentTan && (edgeMu > 0.0 || edgeCohesion > 0.0 || edgeTauMax > 0.0))
        opserr << "WARNING contact -edgeConsistentTan: the non-symmetric edge friction tangent "
                  "needs a non-symmetric solver (system FullGeneral/UmfPack/BandGeneral); "
                  "symmetric solvers will silently corrupt it.\n";
    if (hasMortarOpt && !isMortar) {
        opserr << "WARNING contact -epsN/-epsTie/-epsT/-augTol/-maxAug/-ngp are -mortar "
                  "options; add -mortar or remove them\n";
        return -1;
    }
    if (isMortar) {
        // options of the node-to-segment lane that the mortar lane does not use
        if (consistentNormal) {
            opserr << "WARNING contact -geomtan is a node-to-segment option; the mortar "
                      "tangent does not include the variation of D, M and n\n";
            return -1;
        }
        if (hasCell) {
            opserr << "WARNING contact -cell is a node-to-segment option; mortar pairs every "
                      "slave facet with every master facet\n";
            return -1;
        }
        if (muc > 0.0) {
            opserr << "WARNING contact -visc is a node-to-segment option; it is not available "
                      "with -mortar\n";
            return -1;
        }
        if (isTie && (gapOffset != 0.0 || adjust)) {
            opserr << "WARNING contact -gapOffset/-adjust do not apply to -tie: a tie bonds "
                      "the relative displacement, so the as-meshed gap is already strain-free\n";
            return -1;
        }
        // the refusals of setMortarSmoothing(), checked before the contact is
        // added so that a refused command leaves no definition behind
        if (smoothN > 0.0 || smoothT > 0.0) {
            if (isTie) {
                opserr << "WARNING contact -smoothN/-smoothT do not apply to -tie (a bond has no "
                          "active set and no stick/slip switch to smooth)\n";
                return -1;
            }
            if (augmentMode != LadrunoContactDomain::AUG_NEVER) {
                opserr << "WARNING contact -smoothN/-smoothT are pure-penalty laws; add -augment "
                          "never (the Uzawa update is not defined for the smoothed law)\n";
                return -1;
            }
            if (smoothT > 0.0 && mu <= 0.0 && cohesion <= 0.0 && tauMax <= 0.0) {
                opserr << "WARNING contact -smoothT smooths the stick/slip corner, but this "
                          "contact has no friction (-mu/-cohesion); remove -smoothT\n";
                return -1;
            }
        }
        // the friction options are those of the node-to-segment lane (-mu or
        // the positional mu, -cohesion, -tauMax); a positional kt is the
        // tangential penalty when -epsT is not given
        if (!epsTAuto && epsT <= 0.0 && kt > 0.0)
            epsT = kt;
        int mres = cd->addMortarContact(idata[0], idata[1], idata[2], kn, knAuto, epsN, epsNAuto,
                                        augTol, maxAug, ngp, hasOutward ? outward : 0, isTie,
                                        mu, epsT, epsTAuto, cohesion, tauMax, consistentTan);
        if (mres == 0 && hasControls)
            mres = cd->setMortarContactOptions(idata[0], augmentMode, maxGap, gapOffset,
                                               adjust, adjustTol);
        if (mres == 0 && (smoothN > 0.0 || smoothT > 0.0))
            mres = cd->setMortarSmoothing(idata[0], smoothN, smoothT);
        if (mres == 0 && edgeEdge)
            mres = cd->setMortarEdgeEdge(idata[0], edgeKn, edgeKnAuto, edgeBand, edgeMu, edgeKt,
                                         edgeCohesion, edgeTauMax, edgeConsistentTan, edgeAlm,
                                         edgeAugTol);
        return mres;
    }
    return cd->addContact(idata[0], idata[1], idata[2], kn, kt, mu,
                          hasOutward ? outward : 0, knAuto, cellFrac, consistentTan,
                          consistentNormal, cohesion, tauMax, muc);
}

// contactPlane tag slaveSurfTag nx ny nz px py pz kn <-visc muc>
// A rigid analytical plane through (px, py, pz) with outward normal
// (nx, ny, nz) against the slave node set.
int OPS_LadrunoContactPlane()
{
    if (OPS_GetNumRemainingInputArgs() < 9) {
        opserr << "WARNING want - contactPlane tag slaveSurfTag nx ny nz px py pz kn\n";
        return -1;
    }
    int idata[2], ni = 2;
    if (OPS_GetIntInput(&ni, idata) < 0) {
        opserr << "WARNING contactPlane - could not read tag/slaveSurfTag\n";
        return -1;
    }
    double d[7];
    int nd = 7;
    if (OPS_GetDoubleInput(&nd, d) < 0) {
        opserr << "WARNING contactPlane - could not read nx ny nz px py pz kn\n";
        return -1;
    }
    double n[3]  = { d[0], d[1], d[2] };
    double p0[3] = { d[3], d[4], d[5] };
    double kn = d[6];
    double muc = 0.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (opt != 0 && strcmp(opt, "-visc") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contactPlane -visc - need a coefficient value\n";
                return -1;
            }
            muc = v[0];
        } else if (opt != 0 && strcmp(opt, "-soft") == 0) {
            opserr << "WARNING contactPlane -soft - the soft penalty is sized from the time "
                      "step of an explicit integrator, and no integrator in this build "
                      "reports one; use kn\n";
            return -1;
        } else {
            opserr << "WARNING contactPlane - unexpected token '" << (opt ? opt : "")
                   << "' (expected -visc or end of arguments)\n";
            return -1;
        }
    }
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        opserr << "WARNING contactPlane - define the slave contactSurface first\n";
        return -1;
    }
    return cd->addRigidPlane(idata[0], idata[1], p0, n, kn, muc);
}

// ladrunoContactInfo -> [numContacts, numCommits, numReverts, numMortarContacts]
// (zeros without an engine); numContacts counts the node-to-segment contacts
int OPS_LadrunoContactInfo()
{
    Domain *theDomain = OPS_GetDomain();
    int out[4] = {0, 0, 0, 0};
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0) {
            out[0] = cd->getNumContacts();
            out[1] = cd->getNumCommits();
            out[2] = cd->getNumReverts();
            out[3] = cd->getNumMortarContacts();
        }
    }
    int four = 4;
    if (OPS_SetIntOutput(&four, out, false) < 0)
        return -1;
    return 0;
}

// ladrunoMortarPenetration -> the largest normal penetration max(0, -gbar_I)
// over the mortar slave nodes with a compressive augmented pressure (0 without
// an engine or out of contact). The convergence measure of a held-load
// augmentation loop.
int OPS_LadrunoMortarPenetration()
{
    Domain *theDomain = OPS_GetDomain();
    double pen = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            pen = cd->getMaxMortarPenetration();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &pen, true) < 0)
        return -1;
    return 0;
}

// ladrunoMortarTieResidual -> the largest component of the weighted relative
// displacement |r_I / a_I| over the mortar tie slave nodes (0 without a tie).
// The convergence measure of a held-load augmentation loop on a tie.
int OPS_LadrunoMortarTieResidual()
{
    Domain *theDomain = OPS_GetDomain();
    double res = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            res = cd->getMaxMortarTieResidual();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &res, true) < 0)
        return -1;
    return 0;
}

// ladrunoEdgePenetration -> the largest normal penetration max(0, -gN) over
// the edge-to-edge pairs of -edgeAlm contacts with a compressive augmented
// pressure (0 without such a pair). The convergence measure of a held-load
// augmentation loop on edge-to-edge contact.
int OPS_LadrunoEdgePenetration()
{
    Domain *theDomain = OPS_GetDomain();
    double pen = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            pen = cd->getMaxEdgePenetration();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &pen, true) < 0)
        return -1;
    return 0;
}

// ladrunoBeginAugment / ladrunoEndAugment: open and close a held-load
// augmentation sweep. In between, Domain::commit() updates the contact
// multipliers but records nothing and does not advance the commit tag, so the
// repeated commits of an augmentation loop at a zero load increment add no
// recorder samples. The flag is also cleared by wipe and wipeAnalysis.
int OPS_LadrunoBeginAugment()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain != 0) {
        // While the flag is on, every commit is silent; a second begin without
        // an end means one end was missed and later recorder output is lost.
        if (theDomain->isContactAugmenting())
            opserr << "WARNING ladrunoBeginAugment - an augmentation sweep is already open (a "
                      "ladrunoEndAugment was missed); recorders stay silent until "
                      "ladrunoEndAugment, wipe or wipeAnalysis\n";
        theDomain->setContactAugmenting(true);
    }
    return 0;
}

int OPS_LadrunoEndAugment()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain != 0)
        theDomain->setContactAugmenting(false);
    return 0;
}

// ladrunoContactForce slaveNodeTag -> the normal contact force on a
// node-to-segment slave node, summed over its active master-segment pairs
// (0 without an engine or out of contact). The contact force is assembled by
// the contact adapters, not by elements, so it is not part of nodeReaction.
int OPS_LadrunoContactForce()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING want - ladrunoContactForce slaveNodeTag\n";
        return -1;
    }
    int tag = 0, one = 1;
    if (OPS_GetIntInput(&one, &tag) < 0) {
        opserr << "WARNING ladrunoContactForce - could not read slaveNodeTag\n";
        return -1;
    }
    Domain *theDomain = OPS_GetDomain();
    double f = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            f = cd->getNtsForce(tag);
    }
    if (OPS_SetDoubleOutput(&one, &f, true) < 0)
        return -1;
    return 0;
}
