# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, Jose Antonio Abell. All rights reserved.
"""SplitElastic verification against independent stock quad/brick + dashpots.

Run with the freshly built extension on PYTHONPATH:
  PYTHONPATH=build/Release python3 -m unittest discover -s tests -p test_split_elastic.py -v
NumPy is needed for dense operator comparisons; no output files are kept in the tree.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


SIGNS = np.array([[-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1],
                  [-1,-1,1], [1,-1,1], [1,1,1], [-1,1,1]])


def geometry(d, distorted=True):
    x = (SIGNS[:4 if d == 2 else 8, :d] + 1.) / 2
    if distorted:
        x[1, 0] += .17
        x[2, 1] += .11
        if d == 3:
            x[6] += [.12, -.08, .15]
    return x


def nodal_volumes(x, thickness=1.):
    d = x.shape[1]
    s = SIGNS[:len(x), :d]
    volumes = np.zeros(len(x))
    for xi in s / math.sqrt(3.):
        N = np.prod(1+s*xi, axis=1) / (2**d)
        grad = np.array([[s[i,c]*np.prod([1+s[i,k]*xi[k] for k in range(d) if k != c]) / (2**d)
                          for c in range(d)] for i in range(len(x))])
        volumes += N * np.linalg.det(x.T @ grad) * (thickness if d == 2 else 1.)
    return volumes


def model(d, custom, distorted=True, consistent=False, all_shared=False,
          eta=7., rayleigh=(.17, .003, .005, .002), plane_stress=False):
    """Return canonical (node,dof) ordering, coordinates and independent weights."""
    ops.wipe()
    x = geometry(d, distorted)
    n = len(x)
    shared = [all_shared or abs(xx[0]) < 1.e-12 for xx in x]
    factors = np.linspace(.2, 1.1, n)
    t = .7 if d == 2 else 1.
    weights = nodal_volumes(x, t)
    canonical = []
    ops.model('basic', '-ndm', d, '-ndf', d)
    for i, xx in enumerate(x):
        if not custom or shared[i]:
            ops.node(i+1, *xx)
    if custom:
        ops.model('basic', '-ndm', d, '-ndf', 2*d)
        for i, xx in enumerate(x):
            if not shared[i]:
                ops.node(i+1, *xx)
        args = ['-coupling', eta, '-couplingFactors', *factors]
        if d == 2:
            args += ['-thickness', t, '-planeStress' if plane_stress else '-planeStrain']
        if consistent:
            args += ['-cMass']
        ops.element('SplitElastic', 1, *range(1,n+1), 1234., .23, 2.7, *args)
        for i in range(n):
            canonical += [(i+1, c) for c in range(d if shared[i] else 2*d)]
    else:
        for i, xx in enumerate(x):
            if not shared[i]:
                ops.node(i+1+n, *xx)
        ops.nDMaterial('ElasticIsotropic', 1, 1234./2, .23, 2.7/2)
        for b in range(2):
            ids = [i+1+(n if b and not shared[i] else 0) for i in range(n)]
            if d == 2:
                ops.element('quad', b+1, *ids, t, 'PlaneStress' if plane_stress else 'PlaneStrain', 1)
            else:
                ops.element('stdBrick', b+1, *ids, 1)
        for i in range(n):
            canonical += [(i+1, c) for c in range(d)]
            if not shared[i]:
                canonical += [(i+1+n, c) for c in range(d)]
                ops.uniaxialMaterial('Viscous', i+1, eta*factors[i]*weights[i], 1.)
                ops.element('zeroLength', 10+i, i+1, i+1+n,
                            '-mat', *([i+1]*d), '-dir', *range(1,d+1))
    ops.rayleigh(*rayleigh)
    return canonical, x, weights


def assembled(canonical, m=0., c=0., k=0.):
    ops.wipeAnalysis()
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    # A positive mass term permits the free-body stiffness/damping test to solve.
    ops.integrator('GimmeMCK', 1.+m, c, k)
    ops.algorithm('Linear')
    ops.analysis('Transient')
    assert ops.analyze(1, 0.) == 0
    A = np.array(ops.printA('-ret')).reshape((len(canonical),len(canonical)), order='F')
    permutation = [ops.nodeDOFs(node)[dof] for node, dof in canonical]
    return A[np.ix_(permutation, permutation)]


def transient(d, custom, explicit, steps=600):
    """A source-driven strip/column with shared interface, free/fixed ends.

    In 3-D explicit mode the stock brick reference uses independently assembled
    nodal row-sum mass, since its built-in consistent diagonal is not lumped.
    """
    ops.wipe()
    nx, dx = 8, .125
    section = [(0.,0.),(1.,0.),(1.,1.),(0.,1.)] if d==3 else [(0.,),(1.,)]
    ns = len(section)
    def node(i,j,b=0):
        return 1+i*ns+j+(1000 if b and i>0 else 0)
    def cell(i,b=0):
        if d==2:
            return [node(i,0,b),node(i+1,0,b),node(i+1,1,b),node(i,1,b)]
        return [node(i,0,b),node(i+1,0,b),node(i+1,1,b),node(i,1,b),
                node(i,3,b),node(i+1,3,b),node(i+1,2,b),node(i,2,b)]
    ops.model('basic','-ndm',d,'-ndf',d)
    for j,xx in enumerate(section):
        ops.node(node(0,j),0.,*xx)
    if custom:
        ops.model('basic','-ndm',d,'-ndf',2*d)
    for i in range(1,nx+1):
        for j,xx in enumerate(section):
            ops.node(node(i,j),i*dx,*xx)
            if not custom:
                ops.node(node(i,j,1),i*dx,*xx)
    # Only the transverse y translations are active. Source motion is prescribed
    # at the shared x=0 face; both branch source values are therefore identical.
    for i in range(nx+1):
        for j in range(ns):
            if custom and i>0:
                flags=[1]*(2*d); flags[1]=int(i==nx); flags[d+1]=0
                ops.fix(node(i,j),*flags)
            else:
                flags=[1]*d; flags[1]=int(i==nx)
                ops.fix(node(i,j),*flags)
                if not custom and i>0:
                    flags[1]=0; ops.fix(node(i,j,1),*flags)
    # Constant elastic properties here; graded Rayleigh varies cell by cell.
    E,nu,rho=400.,.2,2.
    nodal_mass = d==3 and explicit and not custom
    if not custom:
        ops.nDMaterial('ElasticIsotropic',1,E/2,nu,0. if nodal_mass else rho/2)
    for i in range(nx):
        factors=[max(0.,1.-(tag-1)//ns/nx) for tag in cell(i)]
        if custom:
            options=['-coupling',10.,'-couplingFactors',*factors]
            if not explicit and d==3: options+=['-cMass']
            ops.element('SplitElastic',i+1,*cell(i),E,nu,rho,*options)
            tags=[i+1]
        else:
            tags=[]
            for b in [0,1]:
                tag=1+i+b*nx; tags.append(tag)
                if d==2: ops.element('quad',tag,*cell(i,b),1.,'PlaneStrain',1)
                else: ops.element('stdBrick',tag,*cell(i,b),1)
        r=i/(nx-1); g=3*r*r-2*r*r*r
        # In nodal-mass reference alphaM acts through nodal Rayleigh damping.
        alpha=0. if nodal_mass else .1*g
        ops.region(i+1,'-ele',*tags,'-rayleigh',alpha,.0003*g,0.,0.)
    for i in range(1,nx+1):
        for j in range(ns):
            volume=dx/ns*(.5 if i==nx else 1.)
            if not custom:
                c=10.*(1-i/nx)*volume
                if c:
                    tag=100+i*ns+j
                    ops.uniaxialMaterial('Viscous',tag,c,1.)
                    ops.element('zeroLength',tag,node(i,j),node(i,j,1),'-mat',tag,'-dir',2)
                if nodal_mass:
                    for b in [0,1]:
                        mass=[rho*volume/2]*d
                        ops.mass(node(i,j,b),*mass)
            # For consistent reference matching, omit alphaM only in 3D explicit
            # on both models; tests still cover all Rayleigh terms in operators.
    if nodal_mass or (d==3 and explicit):
        for i in range(nx):
            r=i/(nx-1);g=3*r*r-2*r*r*r
            tags=[i+1] if custom else [i+1,i+1+nx]
            ops.region(i+1,'-ele',*tags,'-rayleigh',0.,.0003*g,0.,0.)
        # Shared source nodes carry two half contributions in nodal reference.
        if nodal_mass:
            for j in range(ns):ops.mass(node(0,j),*([rho*dx/(2*ns)]*d))
    dt=.0005
    values=[]
    for k in range(steps+2):
        a=math.pi*12.*(k*dt-.08)
        values.append(.01*(1-2*a*a)*math.exp(-a*a))
    ops.timeSeries('Path',1,'-dt',dt,'-values',*values)
    ops.pattern('Plain',1,1)
    for j in range(ns):ops.sp(node(0,j),2,1.)
    ops.constraints('Transformation');ops.numberer('Plain')
    ops.system('Diagonal' if explicit else 'BandGeneral')
    ops.algorithm('Linear')
    ops.integrator('ExplicitDifference') if explicit else ops.integrator('Newmark',.5,.25)
    ops.analysis('Transient')
    histories=[]
    for k in range(steps):
        assert ops.analyze(1,dt)==0
        histories.append([ops.nodeDisp(node(i,0),2) for i in range(nx+1)]+
                         [ops.nodeDisp(node(i,0),d+2) if custom else ops.nodeDisp(node(i,0,1),2)
                          for i in range(1,nx+1)])
    return np.array(histories)


class TestSplitElastic(unittest.TestCase):
    def tearDown(self):
        ops.wipe()

    def test_stock_operator_equivalence(self):
        for d, shared, ps in [(2,False,False),(2,False,True),(2,True,False),
                               (3,False,False),(3,True,False)]:
            with self.subTest(d=d, shared=shared, plane_stress=ps):
                matrices = []
                for custom in [False, True]:
                    order, _, _ = model(d, custom, consistent=d==3, all_shared=shared, plane_stress=ps)
                    M = assembled(order)
                    K = assembled(order,k=1.) - M
                    C = assembled(order,c=1.) - M
                    matrices.append((M,K,C))
                for ref, actual in zip(*matrices):
                    np.testing.assert_allclose(actual, ref, rtol=2.e-12, atol=2.e-12)

    def test_lumping_and_relative_modes(self):
        for d in [2,3]:
            with self.subTest(d=d):
                order, x, weights = model(d,True,consistent=True,rayleigh=(0.,0.,0.,0.))
                Mc = np.array(ops.eleResponse(1,'mass')).reshape(len(order),len(order))
                order, _, _ = model(d,True,consistent=False,rayleigh=(0.,0.,0.,0.))
                M = np.array(ops.eleResponse(1,'mass')).reshape(len(order),len(order))
                C = np.array(ops.eleResponse(1,'damping')).reshape(len(order),len(order))
                np.testing.assert_allclose(M, np.diag(Mc.sum(axis=1)), atol=1.e-14)
                np.testing.assert_allclose(C,C.T,atol=0.)
                self.assertGreaterEqual(np.linalg.eigvalsh(C).min(), -1.e-13)
                common = np.array([.31+.07*c for node,c in order])
                # Both branches carry the same vector at each physical point.
                for i,(node,c) in enumerate(order):
                    common[i]=.31+.07*(c%d)+.013*node
                np.testing.assert_allclose(C@common,0.,atol=1.e-14)
                for (node,c),v in zip(order,common):
                    ops.setNodeVel(node,c+1,v,'-commit')
                self.assertEqual(ops.eleResponse(1,'relativePower')[0],0.)
                v = np.random.default_rng(42).normal(size=len(order))
                for (node,c),value in zip(order,v):
                    ops.setNodeVel(node,c+1,float(value),'-commit')
                power=ops.eleResponse(1,'relativePower')[0]
                self.assertGreater(power,0.)
                self.assertAlmostEqual(power,float(v@C@v),places=12)
                # Region replacement must leave relative damping intact.
                ops.region(1,'-ele',1,'-rayleigh',0.,0.,0.,0.)
                np.testing.assert_allclose(np.array(ops.eleResponse(1,'damping')).reshape(C.shape),C,atol=0.)

    def test_affine_patch_stresses(self):
        for d in [2,3]:
            order,x,_=model(d,True,all_shared=False,eta=0.)
            gradient=np.array([[.01,.03,.02],[-.02,.04,.01],[.03,.01,-.02]])[:d,:d]
            for node in range(1,len(x)+1):
                u=gradient@x[node-1]+.05
                for c in range(d):
                    ops.setNodeDisp(node,c+1,float(u[c]),'-commit')
                    if abs(x[node-1,0])>1.e-12:
                        ops.setNodeDisp(node,d+c+1,float(u[c]),'-commit')
            strain=np.r_[np.diag(gradient),gradient[0,1]+gradient[1,0]]
            if d==3:
                strain=np.r_[strain,gradient[1,2]+gradient[2,1],gradient[2,0]+gradient[0,2]]
            ns=3 if d==2 else 6
            actual=np.array(ops.eleResponse(1,'strains')).reshape(-1,ns)
            np.testing.assert_allclose(actual,np.tile(strain,(len(x),1)),atol=2.e-16)
            mu=1234./(2*(1+.23)); lam=1234.*.23/((1+.23)*(1-2*.23))
            D=np.eye(ns)*mu
            D[:d,:d]=lam+2*mu*np.eye(d)
            stress=D@strain
            np.testing.assert_allclose(np.array(ops.eleResponse(1,'stresses')).reshape(-1,ns),np.tile(stress,(len(x),1)),atol=3.e-13)
            np.testing.assert_allclose(np.array(ops.eleResponse(1,'stresses1'))+.0,
                                       .5*np.array(ops.eleResponse(1,'stresses')),atol=2.e-13)

    def test_database_roundtrip(self):
        for d in [2,3]:
            with self.subTest(d=d), tempfile.TemporaryDirectory() as tmp:
                order,_,_=model(d,True,consistent=d==3)
                for i,(node,c) in enumerate(order):
                    ops.setNodeDisp(node,c+1,math.sin(i)*.01,'-commit')
                    ops.setNodeVel(node,c+1,math.cos(i)*.02,'-commit')
                before={q:ops.eleResponse(1,q) for q in ['stiffness','mass','damping','forces','energy']}
                ops.database('File',str(Path(tmp)/'model'))
                ops.save(1)
                # Clear domain without deleting the active datastore.
                for tag in list(ops.getEleTags()): ops.remove('element',tag)
                for tag in list(ops.getNodeTags()): ops.remove('node',tag)
                ops.restore(1)
                for q,v in before.items():
                    np.testing.assert_allclose(ops.eleResponse(1,q),v,atol=1.e-13,err_msg=q)

    def test_transient_equivalence(self):
        for d in [2,3]:
            for explicit in [False,True]:
                with self.subTest(d=d,explicit=explicit):
                    reference=transient(d,False,explicit)
                    actual=transient(d,True,explicit)
                    np.testing.assert_allclose(actual,reference,rtol=1.e-9,atol=2.e-13)
                    self.assertGreater(float(np.max(np.abs(actual))),.005)

    def test_uniform_excitation_is_common_motion(self):
        for d in [2,3]:
            for consistent in [False,True]:
                with self.subTest(d=d,consistent=consistent):
                    ops.wipe();ops.model('basic','-ndm',d,'-ndf',2*d)
                    x=geometry(d,False)
                    for i,xx in enumerate(x):ops.node(i+1,*xx)
                    options=['-cMass'] if consistent else []
                    ops.element('SplitElastic',1,*range(1,len(x)+1),10.,.2,2.,'-coupling',5.,*options)
                    ops.timeSeries('Constant',1)
                    ops.pattern('UniformExcitation',1,1,'-accel',1)
                    ops.constraints('Plain');ops.numberer('Plain')
                    ops.system('BandGeneral' if consistent else 'Diagonal')
                    ops.algorithm('Linear')
                    ops.integrator('Newmark',.5,.25) if consistent else ops.integrator('ExplicitDifference')
                    ops.analysis('Transient')
                    for step in range(100):
                        self.assertEqual(ops.analyze(1,.001),0)
                        self.assertLess(abs(ops.eleResponse(1,'relativePower')[0]),1.e-25)
                        for i in range(len(x)):
                            v=np.array(ops.nodeVel(i+1))
                            np.testing.assert_allclose(v[:d],v[d:],atol=1.e-13)
                    self.assertLess(ops.nodeVel(1,1),-.09)

    def test_explicit_relative_damping_limit(self):
        # A rigid differential mode has lambda=4*eta/rho. Explicit Euler on
        # velocity is stable only for dt*lambda <= 2, independent of wave CFL.
        dt=.01;rho=2.
        for eta in [60.,120.]:
            with self.subTest(eta=eta):
                ops.wipe();ops.model('basic','-ndm',2,'-ndf',4)
                for i,xx in enumerate(geometry(2,False)):ops.node(i+1,*xx)
                ops.element('SplitElastic',1,1,2,3,4,10.,.2,rho,'-coupling',eta)
                for node in range(1,5):
                    ops.setNodeVel(node,1,1.,'-commit')
                    ops.setNodeVel(node,3,-1.,'-commit')
                ops.constraints('Plain');ops.numberer('Plain');ops.system('Diagonal')
                ops.algorithm('Linear');ops.integrator('ExplicitDifference');ops.analysis('Transient')
                speed=[]
                for step in range(50):
                    self.assertEqual(ops.analyze(1,dt),0)
                    speed.append(abs(ops.nodeVel(1,1)))
                if eta==60.:
                    self.assertLess(max(speed),1.00001)
                    self.assertLess(speed[-1],1.e-8)
                else:
                    self.assertGreater(speed[-1],1.e5)

    def test_invalid_inputs(self):
        def start(ndf=4):
            ops.wipe(); ops.model('basic','-ndm',2,'-ndf',ndf)
            for i,xx in enumerate(geometry(2,False)): ops.node(i+1,*xx)
        cases=[(-1.,.2,1.),(1.,.5,1.),(1.,-.999, -1.),(float('nan'),.2,1.)]
        for values in cases:
            start()
            with self.assertRaises(Exception): ops.element('SplitElastic',1,1,2,3,4,*values)
        for options in [('-coupling',-1.),('-couplingFactors',1.,1.,1.,-1.),
                        ('-lMass','-cMass'),('-planeStress','-planeStrain'),('-unknown',)]:
            start()
            with self.assertRaises(Exception): ops.element('SplitElastic',1,1,2,3,4,10.,.2,1.,*options)
        start(3)
        with self.assertRaises(Exception): ops.element('SplitElastic',1,1,2,3,4,10.,.2,1.)
        start()
        with self.assertRaises(Exception): ops.element('SplitElastic',1,1,4,3,2,10.,.2,1.)
        self.assertEqual(ops.getEleTags(),[])


if __name__=='__main__':
    unittest.main()
