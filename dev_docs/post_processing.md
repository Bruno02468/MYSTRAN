
Exporting Results
=================
There are a few main ways to export data (to verify):

| Location |  Case Control Command  |  Notes       |
| -------- |  --------------------  |  -----       |
| F06/ans  |  DISP=ALL              | Defaults to PRINT         |
| F06/ans  |  DISP(PRINT)=ALL       | Writes to the F06 file (and ANS if DEBUG,200)  |
| OP2      |  DISP(PLOT)=ALL        | Writes to the OP2 file    |
| PUNCH    |  DISP(PUNCH)=ALL       | Writes to the PCH file    |
| CSV      |  DISP(CSV)=ALL         | TODO: Write to a CSV file |
| NEU      |  DISP(NEU)=ALL         | TODO: Writes to a NEU file (currently generated with PARAM,POST,-1) |

The F06 and ANS files are basically the same, so main results support is the same.
However, the ANS file is incomplete, so until they're the same, we'll document it.

OP2 PLOT Support
================

| Command  | Notes   |
| ------   | ------- |
| DISP     |         |
| ACCEL    |         |
| OLOAD    |         |
| SPCFORCE |         |
| MPCFORCE |         |
| STRESS   |         |
| STRAIN   |         |
| FORCE    |         |
| ESE      | N/A     |
| GPFORCE  |         |

Using Sets with Output Requests
===============================

Q: Does Mystran support sets?  For example:
  SET 1 = 1,2,10
  DISP(PRINT) = 1
because it supposedly supports
  DISP(PRINT) = 1,2,10
which begs the question of is the code compatible with sets, when we do?
  SET 1 = 10
  DISP(PRINT) = 1


FAQ
---------
**What is the difference between F06, OP2, NEU, and ANS files?**

**F06** is the main human readable Nastran/Mystran output.  It contains error messages, warnings, etc. which can make things harder to parse.  In general, you'd use this for quick checks, but post-processors like FEMAP / PATRAN / pyNastran can't read it.  An F06 result is requested with a PRINT flag.

The **ANS** is LINK9 (results) only and requires a debug flag to export.  It is simpler to parse, but basically the same as the F06.  An ANS result is requested with a PRINT flag and a ``DEBUG,200,1`` flag.

THe **OP2** is the main binary output file that FEMAP / PATRAN / pyNastran can read.  An OP2 result is requested with a PLOT flag (e.g., ``DISP(PLOT) = ALL``).

The **NEU** is an ASCII file that is supported by FEMAP.  It's not well tested against though.  NEU is selected with the ``PARAM,POST,-1`` flag.

The **PCH** is an even simpler to parse ASCII file that might have some issues writing out.  It's not well tested against though.  A PCH result is requested with a PUNCH flag (e.g., ``DISP(PUNCH) = ALL``).


**Why is my OP2 basically empty (<1 KB)?**

You forgot to request results.  Use DISP(PLOT) = ALL.

**Why is my PCH file not created when I use PUNCH?**

Unlcear.  Seems like a bug.


OP2 Params (TODO)
=================
PARAM,POST,-1
 - default=0 (consistent with Nastran) -> no op2
 - NX/MSC Nastran: activates the OP2 (vs. XDB)
 - Mystran (future option): activates OP2 output (vs. nothing)
 - Mystran current: activates the NEU output (vs. using a flag)

PARAM,OGEOM,YES
 - default=YES (consistent with Nastran)
 - NX/MSC Nastran: activates op2 geometry writing 
 - Mystran current: doesn't exist (NO)
 - Mystran future: YES=activates op2 geometry writing; NO=no geometry


Static Post-Processing Support
==============================

| Result              | F06 | ANS | OP2 | NEU | PCH | CSV |  Notes |
| ------              | --- | --- | --- | --- | --- | --- |  ----- |
| Displacement        | Yes | Yes | Yes | Yes | ??? | No  |        |
| SPC Force           | Yes | Yes | Yes | Yes | ??? | No  |        |
| MPC Force           | Yes | ??? | Yes | ??? | ??? | No  |        |
| Applied Load Vector | Yes | Yes | Yes | ??? | ??? | No  |        |
| Grid Point Weight   | Yes | No  | Yes | ??? | ??? | No  | See Grid Point Weight note |
| Grid Point Force    | Yes | Yes | Yes | ??? | ??? | No  |        |

| Force Result        | F06 | ANS | OP2 | NEU | PCH | CSV |  Notes |
| ------------        | --- | --- | --- | --- | --- | --- |  ----- |
| CELASx              | Yes | Yes | Yes | ??? | ??? | No  |        |
| CROD                | Yes | Yes | Yes | ??? | ??? | No  |        |
| CBUSH               | Yes | Yes | Yes | ??? | ??? | No  |        |
| CBAR                | Yes | Yes | Yes | ??? | ??? | No  |        |
| CBEAM               | N/A | N/A | N/A | ??? | N/A | N/A | card not supported |
| CSHEAR              | Yes | Yes | No  | ??? | ??? | No  | large difference between output structure of Nastran & Mystran |
| CTRIA3 Iso          | Yes | Yes | Yes | ??? | ??? | No  |        |
| TRIA3K/QUAD4K       | Yes | Yes | No  | ??? | ??? | No  |        |
| CQUAD4 Iso-Center   | Yes | Yes | Yes | ??? | ??? | No  |        |
| CQUAD4 Iso-Corner   | No  | No  | No  | ??? | ??? | No  |        |
| CTRIA3 Comp         | No  | No  | No  | ??? | ??? | No  |  [results not calculated](https://github.com/MYSTRANsolver/MYSTRAN/issues/53) |
| CQUAD4 Comp         | No  | No  | No  | ??? | ??? | No  |  [results not calculated](https://github.com/MYSTRANsolver/MYSTRAN/issues/53) |
| Solid               | N/A | N/A | N/A | ??? | N/A | N/A | No outputs (expected) |

| Stress Result       | F06 | ANS | OP2 | NEU | PCH | CSV | Notes                     |
| -------------       | --- | --- | --- | --- | --- | --- | ------------------------- |
| CELASx              | ??? | ??? | ??? | ??? | ??? | No  |                           |
| CROD                | Yes | Yes | Yes | ??? | ??? | No  | no axial/torsion margin   |
| CBUSH               | ??? | ??? | ??? | ??? | ??? | No  |                           |
| CBAR                | Yes | Yes | Yes | ??? | ??? | No  | large difference between output structure of Nastran & Mystran |
| CBEAM               | N/A | N/A | N/A | ??? | N/A | N/A | card not supported        |
| CSHEAR              | Yes | No  | Yes | ??? | ??? | No  |                           |
| CTRIA3 Iso          | Yes | Yes | Yes | ??? | ??? | No  | OP2 plane1/2 results faked; no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)             |
| CQUAD4 Iso-Center   | Yes | Yes | Yes | ??? | ??? | No  | OP2 plane1/2 results faked; no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)             |
| CQUAD4 Iso-Corner   | Yes | Yes | ??? | Yes | ??? | No  |              |
| CTRIA3 Comp         | Yes | Yes | Yes | ??? | ??? | No  | no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)             |
| CQUAD4 Comp         | Yes | Yes | Yes | ??? | ??? | No  | no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)             |
| Solid               | Yes | ??? | Yes | ??? | ??? | No  | Center support only (no corner); No directional vectors; No coordinate system support; no transform support         |


| Strain Result        | F06 | ANS | OP2 | PCH | CSV | Notes                   |
| -------------        | --- | --- | --- | --- | --- | ----------------------- |
| CELASx               | ??? | ??? | ??? | ??? | No  |                         |
| CROD                 | No  | No  | No  | ??? | No  | [1d results not calculated](https://github.com/MYSTRANsolver/MYSTRAN/issues/46) |
| CBUSH                | ??? | ??? | ??? | ??? | No  |                         |
| CBAR                 | No  | No  | No  | No  | No  | [1d results not calculated](https://github.com/MYSTRANsolver/MYSTRAN/issues/46) |
| CBEAM                | N/A | N/A | N/A | N/A | N/A | card not supported      |
| CSHEAR               | Yes | No  | Yes | ??? | No  |                         |
| CTRIA3 Iso           | Yes | Yes | ??? | ??? | No  | OP2 plane1/2 results faked; no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)
| CQUAD4 Iso-Center    | Yes | Yes | Yes | ??? | No  | OP2 plane1/2 results faked; no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)
| CQUAD4 Iso-Corner    | ??? | ??? | Yes | ??? | No  | OP2 plane1/2 results faked; no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)
| CTRIA3 Comp          | Yes | Yes | Yes | ??? | No  | no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)
| CQUAD4 Comp          | Yes | Yes | Yes | ??? | No  | no FIBER/CURV support (FIBER only); no MAXS/MISES support (MISES only)
| Solid                | Yes | ??? | Yes | ??? | No  | Center support only (no corner); No directional vectors; No coordinate system support; no transform support

| Strain Energy Result | F06 | ANS | OP2 | PCH | CSV | Notes |
| -------------------- | --- | --- | --- | --- | --- | ----- |
| CELASx               | No  | No  | No  | No  | No  |       |
| CROD                 | No  | No  | No  | No  | No  |       |
| CBUSH                | No  | No  | No  | No  | No  |       |
| CBAR                 | No  | No  | No  | No  | No  |       |
| CBEAM                | N/A | N/A | N/A | N/A | N/A | card not supported |
| CSHEAR               | No  | No  | No  | No  | No  |       |
| CTRIA3 Iso           | No  | No  | No  | No  | No  |       |
| CQUAD4 Iso-Center    | No  | No  | No  | No  | No  |       |
| CQUAD4 Iso-Corner    | No  | No  | No  | No  | No  |       |
| CTRIA3 Comp          | No  | No  | No  | No  | No  |       |
| CQUAD4 Comp          | No  | No  | No  | No  | No  |       |
| Solid                | No  | No  | No  | No  | No  |       |


Eigen Post-Processing Support
=============================

| Result               | F06 | ANS | OP2 | PCH | CSV | Notes |
| ------               | --- | --- | --- | --- | --- | ----- |
| Grid Point Weight    | Yes | No  | Yes | No  | No  |       |
| Eigenvector          | Yes | ??? | Yes | ??? | No  |       |
| Eigenvalue           | Yes | ??? | Yes | ??? | No  |       |
| MEFFMASS             | ??? | ??? | No  | ??? | No  |       |

Notes
=====

Grid Point Weight Limitation
----------------------------
 ! Generate total mass, first and second moments by summing up mass terms.
 ! XD(i) are components of vector from ref point to a mass point.
 ! At this time, mass units are input units without PARAM WTMASS,
 ! which is what we want for the grid point weight generator.
 ! Later the mass will be converted by multiplying by WTMASS.

 - mass should be (3,1) instead of (1,)
 - CG should be (3,3) instead of (3,1)
 - applies to all outputs

ANS/F06 Limitations
-------------------
 - ANS doesn't support:
   - grid point weight
   - rod strain
   - shear stress/strain

OP2 Limitations
---------------
 - No FIBER/CURV support (FIBER only)
 - No MAXS/MISES support (MISES only)
 - No margin support

PCH Limitations
---------------
1. No GPFORCE support

OP2 specific notes
------------------
1. OP2 follows the roughly MSC 2005r2-2010 output format, which uses the
   XXXXXXXX version code
2. At that time, the output for Mystran compatible results were nearly the
   same between MSC and NX/Simcenter.  Since then, it has diverged
