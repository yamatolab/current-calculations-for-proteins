
StressTensorクラス
~~~~~~~~~~~~~~~~~~
本節では、ストレステンソルクラスについて説明します。
このクラスは、ストレステンソル計算の表式を(ほぼ)そのまま表したメソッドを
使うことが出来ます。
すなわち、 `cal_kinetic` メソッドは速度、体積、速度を与えると、
流れを計算して返します。

.. testsetup:: stress_tensor

   import current
   StressTensor = current.current.StressTensor
   import numpy

例えば、2粒子系の場合を考えます。
各パラメータを以下のように与えると、

.. doctest:: stress_tensor

   >>> vel = numpy.array([[0.11,0.08,-0.02], [0.21,-0.10,0.05]])
   >>> masses  = [1.08, 1.08]
   >>> volumes = [1.2, 1.4]
   >>> scal = StressTensor()
   >>> scal.cal_kinetic(vel, masses, volumes)


同様に、 `cal_bonded` 、 `cal_nonbonded`
によって、いくつかのタイプのポテンシャルについて計算することが出来ます。

先程と同じく2粒子系における計算を以下に示します。

.. doctest:: stress_tensor

   >>> crd = numpy.array([[1.0, 1.5, -0.5], [2.2, -0.1, 0.3]])
   >>> tbfs = numpy.array([1, 1, [0.5, 0.3, 0.2]])

   .. >>> volumes = [1.2, 1.4]
   .. >>> scal = StressTensor()
   .. >>> scal.cal_bonded(crd, tbfs, iszero_tbfs, volumes)




* 最初のセットアップ

  * twobody、topology、setting、グループ辞書を読み込める


* 新しい計算のセットアップ

* 計算を実行する

  * run_atomを実行すると原子毎の流れが計算される。

    * ポテンシャル型別にcurrentを計算することが出来る



* 


currentはポテンシャル型別に計算することが出来る。

.. .. doctest::

..    # for bond, angle, torsion, improper
..    >>> bonded_current = scal.cal_bonded(crd, tbfs, volumes)

..    # for coulomb, vdw
..    >>> nonbonded_current = scal.cal_nonbonded(crd, gen_tbfs, volumes)

..    # for kinetic
..    >>> kinetic_current = scal.cal_kinetic(crd, vel, volumes, masses)



.. .. doctest:: current

..    >>> crd, vel = 
..    >>> cal.run_atom(crd, vel)


以下では、StressTensorクラスの振舞いを記述しています。
