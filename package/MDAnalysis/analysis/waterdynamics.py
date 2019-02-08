# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Water dynamics analysis --- :mod:`MDAnalysis.analysis.waterdynamics`
=======================================================================

:Author: Alejandro Bernardin
:Year: 2014-2015
:Copyright: GNU Public License v3

.. versionadded:: 0.11.0

This module provides functions to analize water dynamics trajectories and water
interactions with other molecules.  The functions in this module are: water
orientational relaxation (WOR) [Yeh1999]_, hydrogen bond lifetimes (HBL)
[Rapaport1983]_, angular distribution (AD) [Grigera1995]_, mean square
displacement (MSD) [Brodka1994]_ and survival probability (SP) [Liu2004]_.

For more information about this type of analysis please refer to
[Araya-Secchi2014]_ (water in a protein cavity) and [Milischuk2011]_ (water in
a nanopore).

.. rubric:: References

.. [Rapaport1983] D.C. Rapaport (1983): Hydrogen bonds in water, Molecular
            Physics: An International Journal at the Interface Between
            Chemistry and Physics, 50:5, 1151-1162.

.. [Yeh1999] Yu-ling Yeh and Chung-Yuan Mou (1999).  Orientational Relaxation
             Dynamics of Liquid Water Studied by Molecular Dynamics Simulation,
             J. Phys. Chem. B 1999, 103, 3699-3705.

.. [Grigera1995] Raul Grigera, Susana G. Kalko and Jorge Fischbarg
                 (1995). Wall-Water Interface.  A Molecular Dynamics Study,
                 Langmuir 1996,12,154-158

.. [Liu2004] Pu Liu, Edward Harder, and B. J. Berne (2004).On the Calculation
             of Diffusion Coefficients in Confined Fluids and Interfaces with
             an Application to the Liquid-Vapor Interface of Water,
             J. Phys. Chem. B 2004, 108, 6595-6602.

.. [Brodka1994] Aleksander Brodka (1994). Diffusion in restricted volume,
                Molecular Physics, 1994, Vol.  82, No. 5, 1075-1078.

.. [Araya-Secchi2014] Araya-Secchi, R., Tomas Perez-Acle, Seung-gu Kang, Tien
                      Huynh, Alejandro Bernardin, Yerko Escalona, Jose-Antonio
                      Garate, Agustin D. Martinez, Isaac E. Garcia, Juan
                      C. Saez, Ruhong Zhou (2014). Characterization of a novel
                      water pocket inside the human Cx26 hemichannel
                      structure. Biophysical journal, 107(3), 599-612.

.. [Milischuk2011] Anatoli A. Milischuk and Branka M. Ladanyi. Structure and
                   dynamics of water confined in silica
                   nanopores. J. Chem. Phys. 135, 174709 (2011); doi:
                   10.1063/1.3657408


Example use of the analysis classes
-----------------------------------

HydrogenBondLifetimes
~~~~~~~~~~~~~~~~~~~~~

Analyzing hydrogen bond lifetimes (HBL) :class:`HydrogenBondLifetimes`, both
continuos and intermittent. In this case we are analyzing how residue 38
interact with a water sphere of radius 6.0 centered on the geometric center of
protein and residue 42. If the hydrogen bond lifetimes are very stable, we can
assume that residue 38 is hydrophilic, on the other hand, if the are very
unstable, we can assume that residue 38 is hydrophobic::

  import MDAnalysis
  from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes as HBL

  u = MDAnalysis.Universe(pdb, trajectory)
  selection1 = "byres name OH2 and sphzone 6.0 protein and resid 42"
  selection2 = "resid 38"
  HBL_analysis = HBL(universe, selection1, selection2, 0, 2000, 30)
  HBL_analysis.run()
  time = 0
  #now we print the data ready to plot. The first two columns are the HBLc vs t
  #plot and the second two columns are the HBLi vs t graph
  for HBLc, HBLi in HBL_analysis.timeseries:
      print("{time} {HBLc} {time} {HBLi}".format(time=time, HBLc=HBLc, HBLi=HBLi))
      time += 1

  #we can also plot our data
  plt.figure(1,figsize=(18, 6))

  #HBL continuos
  plt.subplot(121)
  plt.xlabel('time')
  plt.ylabel('HBLc')
  plt.title('HBL Continuos')
  plt.plot(range(0,time),[column[0] for column in HBL_analysis.timeseries])

  #HBL intermitent
  plt.subplot(122)
  plt.xlabel('time')
  plt.ylabel('HBLi')
  plt.title('HBL Intermitent')
  plt.plot(range(0,time),[column[1] for column in HBL_analysis.timeseries])

  plt.show()

where HBLc is the value for the continuos hydrogen bond lifetimes and HBLi is
the value for the intermittent hydrogen bond lifetime, t0 = 0, tf = 2000 and
dtmax = 30. In this way we create 30 windows timestep (30 values in x
axis). The continuos hydrogen bond lifetimes should decay faster than
intermittent.


WaterOrientationalRelaxation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyzing water orientational relaxation (WOR)
:class:`WaterOrientationalRelaxation`. In this case we are analyzing "how fast"
water molecules are rotating/changing direction. If WOR is very stable we can
assume that water molecules are rotating/changing direction very slow, on the
other hand, if WOR decay very fast, we can assume that water molecules are
rotating/changing direction very fast::

  import MDAnalysis
  from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR

  u = MDAnalysis.Universe(pdb, trajectory)
  selection = "byres name OH2 and sphzone 6.0 protein and resid 42"
  WOR_analysis = WOR(universe, selection, 0, 1000, 20)
  WOR_analysis.run()
  time = 0
  #now we print the data ready to plot. The first two columns are WOR_OH vs t plot,
  #the second two columns are WOR_HH vs t graph and the third two columns are WOR_dip vs t graph
  for WOR_OH, WOR_HH, WOR_dip in WOR_analysis.timeseries:
        print("{time} {WOR_OH} {time} {WOR_HH} {time} {WOR_dip}".format(time=time, WOR_OH=WOR_OH, WOR_HH=WOR_HH,WOR_dip=WOR_dip))
        time += 1

  #now, if we want, we can plot our data
  plt.figure(1,figsize=(18, 6))

  #WOR OH
  plt.subplot(131)
  plt.xlabel('time')
  plt.ylabel('WOR')
  plt.title('WOR OH')
  plt.plot(range(0,time),[column[0] for column in WOR_analysis.timeseries])

  #WOR HH
  plt.subplot(132)
  plt.xlabel('time')
  plt.ylabel('WOR')
  plt.title('WOR HH')
  plt.plot(range(0,time),[column[1] for column in WOR_analysis.timeseries])

  #WOR dip
  plt.subplot(133)
  plt.xlabel('time')
  plt.ylabel('WOR')
  plt.title('WOR dip')
  plt.plot(range(0,time),[column[2] for column in WOR_analysis.timeseries])

  plt.show()

where t0 = 0, tf = 1000 and dtmax = 20. In this way we create 20 windows
timesteps (20 values in the x axis), the first window is created with 1000
timestep average (1000/1), the second window is created with 500 timestep
average(1000/2), the third window is created with 333 timestep average (1000/3)
and so on.

AngularDistribution
~~~~~~~~~~~~~~~~~~~

Analyzing angular distribution (AD) :class:`AngularDistribution` for OH vector,
HH vector and dipole vector. It returns a line histogram with vector
orientation preference. A straight line in the output plot means no
preferential orientation in water molecules. In this case we are analyzing if
water molecules have some orientational preference, in this way we can see if
water molecules are under an electric field or if they are interacting with
something (residue, protein, etc)::

  import MDAnalysis
  from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD

  u = MDAnalysis.Universe(pdb, trajectory)
  selection = "byres name OH2 and sphzone 6.0 (protein and (resid 42 or resid 26) )"
  bins = 30
  AD_analysis = AD(universe,selection,bins)
  AD_analysis.run()
  #now we print data ready to graph. The first two columns are P(cos(theta)) vs cos(theta) for OH vector ,
  #the seconds two columns are P(cos(theta)) vs cos(theta) for HH vector and thirds two columns
  #are P(cos(theta)) vs cos(theta) for dipole vector
  for bin in range(bins):
        print("{AD_analysisOH} {AD_analysisHH} {AD_analysisDip}".format(AD_analysis.graph0=AD_analysis.graph[0][bin], AD_analysis.graph1=AD_analysis.graph[1][bin],AD_analysis.graph2=AD_analysis.graph[2][bin]))

  #and if we want to graph our results
  plt.figure(1,figsize=(18, 6))

  #AD OH
  plt.subplot(131)
  plt.xlabel('cos theta')
  plt.ylabel('P(cos theta)')
  plt.title('PDF cos theta for OH')
  plt.plot([float(column.split()[0]) for column in AD_analysis.graph[0][:-1]],[float(column.split()[1]) for column in AD_analysis.graph[0][:-1]])

  #AD HH
  plt.subplot(132)
  plt.xlabel('cos theta')
  plt.ylabel('P(cos theta)')
  plt.title('PDF cos theta for HH')
  plt.plot([float(column.split()[0]) for column in AD_analysis.graph[1][:-1]],[float(column.split()[1]) for column in AD_analysis.graph[1][:-1]])

  #AD dip
  plt.subplot(133)
  plt.xlabel('cos theta')
  plt.ylabel('P(cos theta)')
  plt.title('PDF cos theta for dipole')
  plt.plot([float(column.split()[0]) for column in AD_analysis.graph[2][:-1]],[float(column.split()[1]) for column in AD_analysis.graph[2][:-1]])

  plt.show()


where `P(cos(theta))` is the angular distribution or angular probabilities.


MeanSquareDisplacement
~~~~~~~~~~~~~~~~~~~~~~

Analyzing mean square displacement (MSD) :class:`MeanSquareDisplacement` for
water molecules. In this case we are analyzing the average distance that water
molecules travels inside protein in XYZ direction (cylindric zone of radius
11[nm], Zmax 4.0[nm] and Zmin -8.0[nm]). A strong rise mean a fast movement of
water molecules, a weak rise mean slow movement of particles::

  import MDAnalysis
  from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

  u = MDAnalysis.Universe(pdb, trajectory)
  selection = "byres name OH2 and cyzone 11.0 4.0 -8.0 protein"
  MSD_analysis = MSD(universe, selection, 0, 1000, 20)
  MSD_analysis.run()
  #now we print data ready to graph. The graph
  #represents MSD vs t
  time = 0
  for msd in MSD_analysis.timeseries:
        print("{time} {msd}".format(time=time, msd=msd))
        time += 1

  #Plot
  plt.xlabel('time')
  plt.ylabel('MSD')
  plt.title('MSD')
  plt.plot(range(0,time),MSD_analysis.timeseries)
  plt.show()


.. _SP-examples:

SurvivalProbability
~~~~~~~~~~~~~~~~~~~

Analyzing survival probability (SP) :class:`SurvivalProbability` for water
molecules. In this case we are analyzing how long water molecules remain in a
sphere of radius 12.3 centered in the geometrical center of resid 42, 26, 34
and 80.  A slow decay of SP means a long permanence time of water molecules in
the zone, on the other hand, a fast decay means a short permanence time::

  import MDAnalysis
  from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP
  import matplotlib.pyplot as plt

  universe = MDAnalysis.Universe(pdb, trajectory)
  selection = "byres name OH2 and sphzone 12.3 (resid 42 or resid 26 or resid 34 or resid 80) "
  sp = SP(universe, selection, verbose=True)
  sp.run(start=0, stop=100, tau_max=20)
  tau_timeseries = sp.tau_timeseries
  sp_timeseries = sp.sp_timeseries

  # print in console
  for tau, sp in zip(tau_timeseries, sp_timeseries):
        print("{time} {sp}".format(time=tau, sp=sp))

  # plot
  plt.xlabel('Time')
  plt.ylabel('SP')
  plt.title('Survival Probability')
  plt.plot(taus, sp_timeseries)
  plt.show()


.. _Output:

Output
------

HydrogenBondLifetimes
~~~~~~~~~~~~~~~~~~~~~

Hydrogen bond lifetimes (HBL) data is returned per window timestep, which is
stored in :attr:`HydrogenBondLifetimes.timeseries` (in all the following
descriptions, # indicates comments that are not part of the output)::

    results = [
        [ # time t0
            <HBL_c>, <HBL_i>
        ],
        [ # time t1
            <HBL_c>, <HBL_i>
        ],
        ...
     ]

WaterOrientationalRelaxation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Water orientational relaxation (WOR) data is returned per window timestep,
which is stored in :attr:`WaterOrientationalRelaxation.timeseries`::

    results = [
        [ # time t0
            <WOR_OH>, <WOR_HH>, <WOR_dip>
        ],
        [ # time t1
            <WOR_OH>, <WOR_HH>, <WOR_dip>
        ],
        ...
     ]

AngularDistribution
~~~~~~~~~~~~~~~~~~~

Angular distribution (AD) data is returned per vector, which is stored in
:attr:`AngularDistribution.graph`. In fact, AngularDistribution returns a
histogram::

    results = [
        [ # OH vector values
          # the values are order in this way: <x_axis  y_axis>
            <cos_theta0 ang_distr0>, <cos_theta1 ang_distr1>, ...
        ],
        [ # HH vector values
            <cos_theta0 ang_distr0>, <cos_theta1 ang_distr1>, ...
        ],
        [ # dip vector values
           <cos_theta0 ang_distr0>, <cos_theta1 ang_distr1>, ...
        ],
     ]

MeanSquareDisplacement
~~~~~~~~~~~~~~~~~~~~~~

Mean Square Displacement (MSD) data is returned in a list, which each element
represents a MSD value in its respective window timestep. Data is stored in
:attr:`MeanSquareDisplacement.timeseries`::

    results = [
         #MSD values orders by window timestep
            <MSD_t0>, <MSD_t1>, ...
     ]

SurvivalProbability
~~~~~~~~~~~~~~~~~~~

Survival Probability (SP) computes two lists: a list of taus (:attr:`SurvivalProbability.tau_timeseries`) and a list of their corresponding mean survival
probabilities (:attr:`SurvivalProbability.sp_timeseries`). Additionally, a list :attr:`SurvivalProbability.sp_timeseries_data` is provided which contains
a list of SPs for each tau, which can be used to compute their distribution, etc.

    results = [ tau1, tau2, ..., tau_n ], [ sp_tau1, sp_tau2, ..., sp_tau_n]

Additionally, for each

Classes
--------

.. autoclass:: HydrogenBondLifetimes
   :members:
   :inherited-members:

.. autoclass:: WaterOrientationalRelaxation
   :members:
   :inherited-members:

.. autoclass:: AngularDistribution
   :members:
   :inherited-members:

.. autoclass:: MeanSquareDisplacement
   :members:
   :inherited-members:

.. autoclass:: SurvivalProbability
   :members:
   :inherited-members:

"""
from __future__ import print_function, division, absolute_import

import warnings

from six.moves import range, zip_longest

import numpy as np
import multiprocessing

import MDAnalysis.analysis.hbonds
from MDAnalysis.lib.log import ProgressMeter


class HydrogenBondLifetimes(object):
    r"""Hydrogen bond lifetime analysis

    This is a autocorrelation function that gives the "Hydrogen Bond Lifetimes"
    (HBL) proposed by D.C. Rapaport [Rapaport1983]_. From this function we can
    obtain the continuous and intermittent behavior of hydrogen bonds in
    time. A fast decay in these parameters indicate a fast change in HBs
    connectivity. A slow decay indicate very stables hydrogen bonds, like in
    ice. The HBL is also know as "Hydrogen Bond Population Relaxation"
    (HBPR). In the continuos case we have:

    .. math::
       C_{HB}^c(\tau) = \frac{\sum_{ij}h_{ij}(t_0)h'_{ij}(t_0+\tau)}{\sum_{ij}h_{ij}(t_0)}

    where :math:`h'_{ij}(t_0+\tau)=1` if there is a H-bond between a pair
    :math:`ij` during time interval :math:`t_0+\tau` (continuos) and
    :math:`h'_{ij}(t_0+\tau)=0` otherwise. In the intermittent case we have:

    .. math::
       C_{HB}^i(\tau) = \frac{\sum_{ij}h_{ij}(t_0)h_{ij}(t_0+\tau)}{\sum_{ij}h_{ij}(t_0)}

    where :math:`h_{ij}(t_0+\tau)=1` if there is a H-bond between a pair
    :math:`ij` at time :math:`t_0+\tau` (intermittent) and
    :math:`h_{ij}(t_0+\tau)=0` otherwise.


    Parameters
    ----------
    universe : Universe
      Universe object
    selection1 : str
      Selection string for first selection [‘byres name OH2’].
      It could be any selection available in MDAnalysis, not just water.
    selection2 : str
      Selection string to analize its HBL against selection1
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.
    nproc : int
      Number of processors to use, by default is 1.


    .. versionadded:: 0.11.0
    """

    def __init__(self, universe, selection1, selection2, t0, tf, dtmax,
                 nproc=1):
        self.universe = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.t0 = t0
        self.tf = tf - 1
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _getC_i(self, HBP, t0, t):
        """
        This function give the intermitent Hydrogen Bond Lifetime
        C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
        """
        C_i = 0
        for i in range(len(HBP[t0])):
            for j in range(len(HBP[t])):
                if (HBP[t0][i][0] == HBP[t][j][0] and HBP[t0][i][1] == HBP[t][j][1]):
                    C_i += 1
                    break
        if len(HBP[t0]) == 0:
            return 0.0
        else:
            return float(C_i) / len(HBP[t0])

    def _getC_c(self, HBP, t0, t):
        """
        This function give the continous Hydrogen Bond Lifetime
        C_c = <h(t0)h'(t)>/<h(t0)> between t0 and t
        """
        C_c = 0
        dt = 1
        begt0 = t0
        HBP_cp = HBP
        HBP_t0 = HBP[t0]
        newHBP = []
        if t0 == t:
            return 1.0
        while t0 + dt <= t:
            for i in range(len(HBP_t0)):
                for j in range(len(HBP_cp[t0 + dt])):
                    if (HBP_t0[i][0] == HBP_cp[t0 + dt][j][0] and
                        HBP_t0[i][1] == HBP_cp[t0 + dt][j][1]):
                        newHBP.append(HBP_t0[i])
                        break
            C_c = len(newHBP)
            t0 += dt
            HBP_t0 = newHBP
            newHBP = []
        if len(HBP[begt0]) == 0:
            return 0
        else:
            return C_c / float(len(HBP[begt0]))

    def _intervC_c(self, HBP, t0, tf, dt):
        """
        This function gets all the data for the h(t0)h(t0+dt)', where
        t0 = 1,2,3,...,tf. This function give us one point of the final plot
        HBL vs t
        """
        a = 0
        count = 0
        for i in range(len(HBP)):
            if (t0 + dt <= tf):
                if t0 == t0 + dt:
                    b = self._getC_c(HBP, t0, t0)
                    break
                b = self._getC_c(HBP, t0, t0 + dt)
                t0 += dt
                a += b
                count += 1
        if count == 0:
            return 1.0
        return a / count

    def _intervC_i(self, HBP, t0, tf, dt):
        """
        This function gets all the data for the h(t0)h(t0+dt), where
        t0 = 1,2,3,...,tf. This function give us a point of the final plot
        HBL vs t
        """
        a = 0
        count = 0
        for i in range(len(HBP)):
            if (t0 + dt <= tf):
                b = self._getC_i(HBP, t0, t0 + dt)
                t0 += dt
                a += b
                count += 1
        return a / count

    def _finalGraphGetC_i(self, HBP, t0, tf, maxdt):
        """
        This function gets the final data of the C_i graph.
        """
        output = []
        for dt in range(maxdt):
            a = self._intervC_i(HBP, t0, tf, dt)
            output.append(a)
        return output

    def _finalGraphGetC_c(self, HBP, t0, tf, maxdt):
        """
        This function gets the final data of the C_c graph.
        """
        output = []
        for dt in range(maxdt):
            a = self._intervC_c(HBP, t0, tf, dt)
            output.append(a)
        return output

    def _getGraphics(self, HBP, t0, tf, maxdt):
        """
        Function that join all the results into a plot.
        """
        a = []
        cont = self._finalGraphGetC_c(HBP, t0, tf, maxdt)
        inte = self._finalGraphGetC_i(HBP, t0, tf, maxdt)
        for i in range(len(cont)):
            fix = [cont[i], inte[i]]
            a.append(fix)
        return a

    def _HBA(self, ts, conn, universe, selAtom1, selAtom2,
             verbose=False):
        """
        Main function for calculate C_i and C_c in parallel.
        """
        finalGetResidue1 = selAtom1
        finalGetResidue2 = selAtom2
        frame = ts.frame
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(universe,
                                                            finalGetResidue1,
                                                            finalGetResidue2,
                                                            distance=3.5,
                                                            angle=120.0,
                                                            start=frame - 1,
                                                            stop=frame)
        while True:
            try:
                h.run(verbose=verbose)
                break
            except:
                print("error")
                print("trying again")
                sys.stdout.flush()
        sys.stdout.flush()
        conn.send(h.timeseries[0])
        conn.close()

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""
        h_list = []
        i = 0
        if (self.nproc > 1):
            while i < len(self.universe.trajectory):
                jobs = []
                k = i
                for j in range(self.nproc):
                        # start
                    print("ts=", i + 1)
                    if i >= len(self.universe.trajectory):
                        break
                    conn_parent, conn_child = multiprocessing.Pipe(False)
                    while True:
                        try:
                            # new thread
                            jobs.append(
                                (multiprocessing.Process(
                                    target=self._HBA,
                                    args=(self.universe.trajectory[i],
                                          conn_child, self.universe,
                                          self.selection1, self.selection2,)),
                                 conn_parent))
                            break
                        except:
                            print("error in jobs.append")
                    jobs[j][0].start()
                    i = i + 1

                for j in range(self.nproc):
                    if k >= len(self.universe.trajectory):
                        break
                    rec01 = jobs[j][1]
                    received = rec01.recv()
                    h_list.append(received)
                    jobs[j][0].join()
                    k += 1
            self.timeseries = self._getGraphics(
                h_list, 0, self.tf - 1, self.dtmax)
        else:
            h_list = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe,
                                                                     self.selection1,
                                                                     self.selection2,
                                                                     distance=3.5,
                                                                     angle=120.0)
            h_list.run(**kwargs)
            self.timeseries = self._getGraphics(
                h_list.timeseries, self.t0, self.tf, self.dtmax)


class WaterOrientationalRelaxation(object):
    r"""Water orientation relaxation analysis

    Function to evaluate the Water Orientational Relaxation proposed by Yu-ling
    Yeh and Chung-Yuan Mou [Yeh1999_]. WaterOrientationalRelaxation indicates
    "how fast" water molecules are rotating or changing direction. This is a
    time correlation function given by:

    .. math::
        C_{\hat u}(\tau)=\langle \mathit{P}_2[\mathbf{\hat{u}}(t_0)\cdot\mathbf{\hat{u}}(t_0+\tau)]\rangle

    where :math:`P_2=(3x^2-1)/2` is the second-order Legendre polynomial and :math:`\hat{u}` is
    a unit vector along HH, OH or dipole vector.


    Parameters
    ----------
    universe : Universe
      Universe object
    selection : str
      Selection string for water [‘byres name OH2’].
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.


    .. versionadded:: 0.11.0

    """

    def __init__(self, universe, selection, t0, tf, dtmax, nproc=1):
        self.universe = universe
        self.selection = selection
        self.t0 = t0
        self.tf = tf
        # make sure range specified is appropriate
        self.universe.trajectory[tf-1]
        self.total_frames = self.tf - self.t0
        
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _get_water_vectors(self, frame, selection):
        """
        Calculates the vectors along each water molecule in selection
        in the given frame.
        Returns: A 1x3 np array. Think of it as a list of all the 3D vectors.
        """
        self.universe.trajectory[frame]

        # The selection has atoms arranged as [O H1 H2 O H1 H2...]
        oxygens = np.array([atom.position for atom in selection[0::3]])
        h1s = np.array([atom.position for atom in selection[1::3]])
        h2s = np.array([atom.position for atom in selection[2::3]])

        OH_vecs = h1s - oxygens
        HH_vecs = h1s - h2s
        DIP_vecs = (h1s + h2s) * 0.5 - oxygens

        return (OH_vecs, HH_vecs, DIP_vecs)

    def _get_one_delta_point(self, selection, t0, t1):
        """
        Calculates the dot product of each UNIT water vector at t0 and t1,
        then takes the Lg polynomial of that result, and avegares over all the
        water molecules in the current frames being compared.

        In calculating the normalized dot product, it is faster to calculate
        U*V / sqrt(||U||^2 x ||V||^2) than it is to calculate U/||U|| * V/||V||.
        Note that ||U||^2 = U*U.

        einsum('ij,ij->i', A, B) will multiply corresponding entries in both
          arrays before summing down the columns to yield a 1D array. When
          applied to lists of vectors as we have here, the result is the list of
          dot products of corresponding vectors: [A[0]*B[0], A[1]*B[1], ...].
          Einsum was faster in this scenario than calling np.linalg.norm on each
          entry.
        """
        OHi, HHi, DIi = self._get_water_vectors(t0, selection)
        OHf, HHf, DIf = self._get_water_vectors(t1, selection)

        dot = 'ij,ij->i'

        # Compute the product of the norms
        OH_norms = np.sqrt(np.einsum(dot, OHi, OHi) * np.einsum(dot, OHf, OHf))
        HH_norms = np.sqrt(np.einsum(dot, HHi, HHi) * np.einsum(dot, HHf, HHf))
        DI_norms = np.sqrt(np.einsum(dot, DIi, DIi) * np.einsum(dot, DIf, DIf))

        # Compute the dot products and normalize
        OH_dots = np.einsum(dot, OHi, OHf) / OH_norms
        HH_dots = np.einsum(dot, HHi, HHf) / HH_norms
        DI_dots = np.einsum(dot, DIi, DIf) / DI_norms

        OH_vals = self.lg2(OH_dots)
        HH_vals = self.lg2(HH_dots)
        DI_vals = self.lg2(DI_dots)

        return np.array([np.mean(OH_vals), np.mean(HH_vals), np.mean(DI_vals)])


    def _get_mean_one_point(self, selection1, dt):
        """
        This function gets one point of the plot C_vec vs t. It uses the
        _get_one_delta_point() function to calculate the average.

        """
        repInd = self._repeatedIndex(selection1, dt)
        n = 0.0
        delta_sums = np.zeros(3)

        for index, t_cur in enumerate(range(self.t0, self.tf, dt)):
            t_next = t_cur + dt
            if t_next >= self.tf:
                break
            delta_sums += self._get_one_delta_point(repInd[index], t_cur, t_next)
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return delta_sums / n if n > 0 else (0, 0, 0)

    def _repeatedIndex(self, selection, dt):
        """
        Indicates the comparation between all the t+dt.
        The results is a list of list with all the repeated index per frame
        (or time).
        Ex: dt=1, so compare frames (1,2),(2,3),(3,4)...
        Ex: dt=2, so compare frames (1,3),(3,5),(5,7)...
        Ex: dt=3, so compare frames (1,4),(4,7),(7,10)...
        """
        rep = []
        for i in range(0, self.total_frames, dt):
            j = i + dt
            if j >= self.total_frames:
                break
            rep.append(self._sameMolecTandDT(selection, i, j))
        return rep

    def _sameMolecTandDT(self, selection, t0d, tf):
        """
        Compare the molecules in the t0d selection and the t0d+dt selection and
        select only the particles that are repeated in both frame. This is to
        consider only the molecules that remains in the selection after the dt
        time has elapsed.
        The result is a list with the indexs of the atoms.
        """
        a = set(selection[t0d])
        b = set(selection[tf])
        sort = sorted(list(a.intersection(b)))
        return sort

    def _selection_serial(self, universe, selection_str):
        selection = []
        pm = ProgressMeter(self.total_frames,
                           interval=10, verbose=True)
        for ts in universe.trajectory[self.t0:self.tf]:
            selection.append(universe.select_atoms(selection_str))
            pm.echo(ts.frame - self.t0)
        return selection

    @staticmethod
    def lg2(x):
        """Second Legendre polynomial"""
        return (3*x*x - 1)/2

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""

        # All the selection to an array, this way is faster than selecting
        # later.
        if self.nproc == 1:
            selection_out = self._selection_serial(
                self.universe, self.selection)
        else:
            # selection_out = self._selection_parallel(self.universe,
            # self.selection, self.nproc)
            # parallel selection to be implemented
            selection_out = self._selection_serial(
                self.universe, self.selection)
        self.timeseries = []
        for dt in list(range(1, self.dtmax + 1)):
            output = self._get_mean_one_point(selection_out, dt)
            self.timeseries.append(output)


class AngularDistribution(object):
    r"""Angular distribution function analysis

    The angular distribution function (AD) is defined as the distribution
    probability of the cosine of the :math:`\theta` angle formed by the OH
    vector, HH vector or dipolar vector of water molecules and a vector
    :math:`\hat n` parallel to chosen axis (z is the default value). The cosine
    is define as :math:`\cos \theta = \hat u \cdot \hat n`, where :math:`\hat
    u` is OH, HH or dipole vector.  It creates a histogram and returns a list
    of lists, see Output_. The AD is also know as Angular Probability (AP).


    Parameters
    ----------
    universe : Universe
        Universe object
    selection : str
        Selection string to evaluate its angular distribution ['byres name OH2']
    bins : int (optional)
        Number of bins to create the histogram by means of :func:`numpy.histogram`
    axis : {'x', 'y', 'z'} (optional)
        Axis to create angle with the vector (HH, OH or dipole) and calculate
        cosine theta ['z'].


    .. versionadded:: 0.11.0
    """

    def __init__(self, universe, selection_str, bins=40, nproc=1, axis="z"):
        self.universe = universe
        self.selection_str = selection_str
        self.bins = bins
        self.nproc = nproc
        self.axis = axis
        self.graph = None

    def _getCosTheta(self, universe, selection, axis):
        valOH = []
        valHH = []
        valdip = []

        i = 0
        while i <= (len(selection) - 1):
            universe.trajectory[i]
            line = selection[i].positions

            Ot0 = line[::3]
            H1t0 = line[1::3]
            H2t0 = line[2::3]

            OHVector0 = H1t0 - Ot0
            HHVector0 = H1t0 - H2t0
            dipVector0 = (H1t0 + H2t0) * 0.5 - Ot0

            unitOHVector0 = OHVector0 / \
                np.linalg.norm(OHVector0, axis=1)[:, None]
            unitHHVector0 = HHVector0 / \
                np.linalg.norm(HHVector0, axis=1)[:, None]
            unitdipVector0 = dipVector0 / \
                np.linalg.norm(dipVector0, axis=1)[:, None]

            j = 0
            while j < len(line) / 3:
                if axis == "z":
                    valOH.append(unitOHVector0[j][2])
                    valHH.append(unitHHVector0[j][2])
                    valdip.append(unitdipVector0[j][2])

                elif axis == "x":
                    valOH.append(unitOHVector0[j][0])
                    valHH.append(unitHHVector0[j][0])
                    valdip.append(unitdipVector0[j][0])

                elif axis == "y":
                    valOH.append(unitOHVector0[j][1])
                    valHH.append(unitHHVector0[j][1])
                    valdip.append(unitdipVector0[j][1])

                j += 1
            i += 1
        return (valOH, valHH, valdip)

    def _getHistogram(self, universe, selection, bins, axis):
        """
        This function gets a normalized histogram of the cos(theta) values. It
        return a list of list.
        """
        a = self._getCosTheta(universe, selection, axis)
        cosThetaOH = a[0]
        cosThetaHH = a[1]
        cosThetadip = a[2]
        lencosThetaOH = len(cosThetaOH)
        lencosThetaHH = len(cosThetaHH)
        lencosThetadip = len(cosThetadip)
        histInterval = bins
        histcosThetaOH = np.histogram(cosThetaOH, histInterval, density=True)
        histcosThetaHH = np.histogram(cosThetaHH, histInterval, density=True)
        histcosThetadip = np.histogram(cosThetadip, histInterval, density=True)

        return (histcosThetaOH, histcosThetaHH, histcosThetadip)

    def _hist2column(self, aList):
        """
        This function transform from the histogram format
        to a column format.
        """
        a = []
        for x in zip_longest(*aList, fillvalue="."):
            a.append(" ".join(str(i) for i in x))
        return a

    def run(self, **kwargs):
        """Function to evaluate the angular distribution of cos(theta)"""

        if self.nproc == 1:
            selection = self._selection_serial(
                self.universe, self.selection_str)
        else:
            # not implemented yet
            # selection = self._selection_parallel(self.universe,
            # self.selection_str,self.nproc)
            selection = self._selection_serial(
                self.universe, self.selection_str)

        self.graph = []
        output = self._getHistogram(
            self.universe, selection, self.bins, self.axis)
        # this is to format the exit of the file
        # maybe this output could be improved
        listOH = [list(output[0][1]), list(output[0][0])]
        listHH = [list(output[1][1]), list(output[1][0])]
        listdip = [list(output[2][1]), list(output[2][0])]

        self.graph.append(self._hist2column(listOH))
        self.graph.append(self._hist2column(listHH))
        self.graph.append(self._hist2column(listdip))

    def _selection_serial(self, universe, selection_str):
        selection = []
        pm = ProgressMeter(universe.trajectory.n_frames,
                           interval=10, verbose=True)
        for ts in universe.trajectory:
            selection.append(universe.select_atoms(selection_str))
            pm.echo(ts.frame)
        return selection


class MeanSquareDisplacement(object):
    r"""Mean square displacement analysis

    Function to evaluate the Mean Square Displacement (MSD_). The MSD gives the
    average distance that particles travels. The MSD is given by:

    .. math::
        \langle\Delta r(t)^2\rangle = 2nDt

    where :math:`r(t)` is the position of particle in time :math:`t`,
    :math:`\Delta r(t)` is the displacement after time lag :math:`t`,
    :math:`n` is the dimensionality, in this case :math:`n=3`,
    :math:`D` is the diffusion coefficient and :math:`t` is the time.

    .. _MSD: http://en.wikipedia.org/wiki/Mean_squared_displacement


    Parameters
    ----------
    universe : Universe
      Universe object
    selection : str
      Selection string for water [‘byres name OH2’].
    t0 : int
      frame  where analysis begins
    tf : int
      frame where analysis ends
    dtmax : int
      Maximum dt size, `dtmax` < `tf` or it will crash.


    .. versionadded:: 0.11.0
    """

    def __init__(self, universe, selection, t0, tf, dtmax, nproc=1):
        self.universe = universe
        self.selection = selection
        self.t0 = t0
        self.tf = tf
        self.dtmax = dtmax
        self.nproc = nproc
        self.timeseries = None

    def _repeatedIndex(self, selection, dt, totalFrames):
        """
        Indicate the comparation between all the t+dt.
        The results is a list of list with all the repeated index per frame
        (or time).

        - Ex: dt=1, so compare frames (1,2),(2,3),(3,4)...
        - Ex: dt=2, so compare frames (1,3),(3,5),(5,7)...
        - Ex: dt=3, so compare frames (1,4),(4,7),(7,10)...
        """
        rep = []
        for i in range(int(round((totalFrames - 1) / float(dt)))):
            if (dt * i + dt < totalFrames):
                rep.append(self._sameMolecTandDT(
                    selection, dt * i, (dt * i) + dt))
        return rep

    def _getOneDeltaPoint(self, universe, repInd, i, t0, dt):
        """
        Gives one point to calculate the mean and gets one point of the plot
        C_vect vs t.

        - Ex: t0=1 and dt=1 so calculate the t0-dt=1-2 interval.
        - Ex: t0=5 and dt=3 so calcultate the t0-dt=5-8 interva

        i = come from getMeanOnePoint (named j) (int)
        """
        valO = 0
        n = 0
        for j in range(len(repInd[i]) // 3):
            begj = 3 * j
            universe.trajectory[t0]
            # Plus zero is to avoid 0to be equal to 0tp
            Ot0 = repInd[i][begj].position + 0

            universe.trajectory[t0 + dt]
            # Plus zero is to avoid 0to be equal to 0tp
            Otp = repInd[i][begj].position + 0

            # position oxygen
            OVector = Ot0 - Otp
            # here it is the difference with
            # waterdynamics.WaterOrientationalRelaxation
            valO += np.dot(OVector, OVector)
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return valO/n if n > 0 else 0

    def _getMeanOnePoint(self, universe, selection1, selection_str, dt,
                         totalFrames):
        """
        This function gets one point of the plot C_vec vs t. It's uses the
        _getOneDeltaPoint() function to calculate the average.

        """
        repInd = self._repeatedIndex(selection1, dt, totalFrames)
        sumsdt = 0
        n = 0.0
        sumDeltaO = 0.0
        valOList = []

        for j in range(totalFrames // dt - 1):
            a = self._getOneDeltaPoint(universe, repInd, j, sumsdt, dt)
            sumDeltaO += a
            valOList.append(a)
            sumsdt += dt
            n += 1

        # if no water molecules remain in selection, there is nothing to get
        # the mean, so n = 0.
        return sumDeltaO/n if n > 0 else 0

    def _sameMolecTandDT(self, selection, t0d, tf):
        """
        Compare the molecules in the t0d selection and the t0d+dt selection and
        select only the particles that are repeated in both frame. This is to
        consider only the molecules that remains in the selection after the dt
        time has elapsed. The result is a list with the indexs of the atoms.
        """
        a = set(selection[t0d])
        b = set(selection[tf])
        sort = sorted(list(a.intersection(b)))
        return sort

    def _selection_serial(self, universe, selection_str):
        selection = []
        pm = ProgressMeter(universe.trajectory.n_frames,
                           interval=10, verbose=True)
        for ts in universe.trajectory:
            selection.append(universe.select_atoms(selection_str))
            pm.echo(ts.frame)
        return selection

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries"""

        # All the selection to an array, this way is faster than selecting
        # later.
        if self.nproc == 1:
            selection_out = self._selection_serial(
                self.universe, self.selection)
        else:
            # parallel not yet implemented
            # selection = selection_parallel(universe, selection_str, nproc)
            selection_out = self._selection_serial(
                self.universe, self.selection)
        self.timeseries = []
        for dt in list(range(1, self.dtmax + 1)):
            output = self._getMeanOnePoint(
                self.universe, selection_out, self.selection, dt, self.tf)
            self.timeseries.append(output)


class SurvivalProbability(object):
    r"""Survival probability analysis

    Function to evaluate the Survival Probability (SP). The SP gives the
    probability for a group of particles to remain in certain region. The SP is
    given by:

    .. math::
        P(\tau) = \frac1T \sum_{t=1}^T \frac{N(t,t+\tau)}{N(t)}

    where :math:`T` is the maximum time of simulation, :math:`\tau` is the
    timestep, :math:`N(t)` the number of particles at time t, and
    :math:`N(t, t+\tau)` is the number of particles at every frame from t to `\tau`.


    Parameters
    ----------
    universe : Universe
      Universe object
    selection : str
      Selection string; any selection is allowed. With this selection you
      define the region/zone where to analyze, e.g.: "resname SOL and around 5 (resname LIPID)"
      and "resname ION and around 10 (resid 20)" (see `SP-examples`_ )
    verbose : Boolean
      If True, prints progress and comments to the console.


    .. versionadded:: 0.11.0

    """

    def __init__(self, universe, selection, t0=None, tf=None, dtmax=None, verbose=False):
        self.universe = universe
        self.selection = selection
        self.verbose = verbose

        # backward compatibility
        self.start = self.stop = self.tau_max = None
        if t0 is not None:
            self.start = t0
            warnings.warn("t0 is deprecated, use run(start=t0) instead", category=DeprecationWarning)

        if tf is not None:
            self.stop = tf
            warnings.warn("tf is deprecated, use run(stop=tf) instead", category=DeprecationWarning)

        if dtmax is not None:
            self.tau_max = dtmax
            warnings.warn("dtmax is deprecated, use run(tau_max=dtmax) instead", category=DeprecationWarning)

    def print(self, verbose, *args):
        if self.verbose:
            print(args)
        elif verbose:
            print(args)

    def run(self, tau_max=20, start=0, stop=None, step=1, verbose=False):
        """
        Computes and returns the survival probability timeseries

        Parameters
        ----------
        start : int
            Zero-based index of the first frame to be analysed
        stop : int
            Zero-based index of the last frame to be analysed (inclusive)
        step : int
            Jump every `step`'th frame
        tau_max : int
            Survival probability is calculated for the range :math:`1 <= \tau <= tau_max`
        verbose : Boolean
            Overwrite the constructor's verbosity

        Returns
        -------
        tau_timeseries : list
            tau from 1 to tau_max. Saved in the field tau_timeseries.
        sp_timeseries : list
            survival probability for each value of `tau`. Saved in the field sp_timeseries.
        """

        # backward compatibility (and priority)
        start = self.start if self.start is not None else start
        stop = self.stop if self.stop is not None else stop
        tau_max = self.tau_max if self.tau_max is not None else tau_max

        # sanity checks
        if stop is not None and stop >= len(self.universe.trajectory):
            raise ValueError("\"stop\" must be smaller than the number of frames in the trajectory.")

        if stop is None:
            stop = len(self.universe.trajectory)
        else:
            stop = stop + 1

        if tau_max > (stop - start):
            raise ValueError("Too few frames selected for given tau_max.")

        # load all frames to an array of sets
        selected_ids = []
        for ts in self.universe.trajectory[start:stop]:
            self.print(verbose, "Loading frame:", ts)
            selected_ids.append(set(self.universe.select_atoms(self.selection).ids))

        tau_timeseries = np.arange(1, tau_max + 1)
        sp_timeseries_data = [[] for _ in range(tau_max)]

        for t in range(0, len(selected_ids), step):
            Nt = len(selected_ids[t])

            if Nt == 0:
                self.print(verbose,
                           "At frame {} the selection did not find any molecule. Moving on to the next frame".format(t))
                continue

            for tau in tau_timeseries:
                if t + tau >= len(selected_ids):
                    break

                # ids that survive from t to t + tau and at every frame in between
                Ntau = len(set.intersection(*selected_ids[t:t + tau + 1]))
                sp_timeseries_data[tau - 1].append(Ntau / float(Nt))

        # user can investigate the distribution and sample size
        self.sp_timeseries_data = sp_timeseries_data

        self.tau_timeseries = tau_timeseries
        self.sp_timeseries = [np.mean(sp) for sp in sp_timeseries_data]
        return self
