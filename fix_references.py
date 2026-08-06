"""G-4: replace every bibliography entry with one verified against the publisher record.

Each entry below was checked on 6 Aug 2026 against the publisher's own record or,
for the datasets, the repository landing page. Corrections found:

  richardson2015 - DOI was missing
  lin2013        - volume, issue, page range and full author list were missing
  perez2012      - TITLE WAS WRONG: the published title contains "cylindrical",
                   which the source archive's readme omits. Page range, paper
                   number, venue detail and DOI were all missing.
  forgez2010     - DOI was missing
  catenaro2021   - the companion Data in Brief article was not cited. The article
                   references Mendeley kxsbr4x3j2 specifically, and lists exactly
                   the eight LFP C-rates used here, so it documents the data
                   actually analysed.
  bernardi1985   - already correct
  saha2007       - repository entry, no DOI issued
"""
import io
import re

VERIFIED = {
    "richardson2015": r"""R.~R. Richardson and D.~A. Howey, ``Sensorless battery internal temperature
estimation using a Kalman filter with impedance measurement,'' \emph{IEEE Trans.
Sustain. Energy}, vol.~6, no.~4, pp.~1190--1199, Oct. 2015,
doi: 10.1109/TSTE.2015.2420375.""",

    "catenaro2021": r"""E.~Catenaro and S.~Onori, ``Experimental data of three lithium-ion batteries
under galvanostatic discharge tests at different C-rates and operating
temperatures,'' Mendeley Data, 2021, doi: 10.17632/kxsbr4x3j2.1. Documented in
E.~Catenaro and S.~Onori, ``Experimental data of lithium-ion batteries under
galvanostatic discharge tests at different rates and temperatures of
operation,'' \emph{Data in Brief}, vol.~35, art.~106894, 2021,
doi: 10.1016/j.dib.2021.106894.""",

    "saha2007": r"""B.~Saha and K.~Goebel, ``Battery data set,'' NASA Ames Prognostics Data
Repository, NASA Ames Research Center, Moffett Field, CA, USA, 2007.""",

    "lin2013": r"""X.~Lin, H.~E. Perez, J.~B. Siegel, A.~G. Stefanopoulou, Y.~Li, R.~D. Anderson,
Y.~Ding, and M.~P. Castanier, ``Online parameterization of lumped thermal
dynamics in cylindrical lithium ion batteries for core temperature estimation
and health monitoring,'' \emph{IEEE Trans. Control Syst. Technol.}, vol.~21,
no.~5, pp.~1745--1755, Sep. 2013, doi: 10.1109/TCST.2012.2217143.""",

    "perez2012": r"""H.~E. Perez, J.~B. Siegel, X.~Lin, A.~G. Stefanopoulou, Y.~Ding, and M.~P.
Castanier, ``Parameterization and validation of an integrated electro-thermal
cylindrical LFP battery model,'' in \emph{Proc. ASME 5th Annu. Dyn. Syst.
Control Conf. joint JSME 11th Motion Vib. Conf.}, Fort Lauderdale, FL, USA,
Oct. 2012, pp.~41--50, paper DSCC2012-MOVIC2012-8782,
doi: 10.1115/DSCC2012-MOVIC2012-8782.""",

    "bernardi1985": r"""D.~Bernardi, E.~Pawlikowski, and J.~Newman, ``A general energy balance for
battery systems,'' \emph{J. Electrochem. Soc.}, vol.~132, no.~1, pp.~5--12,
1985, doi: 10.1149/1.2113792.""",

    "forgez2010": r"""C.~Forgez, D.~V. Do, G.~Friedrich, M.~Morcrette, and C.~Delacourt, ``Thermal
modeling of a cylindrical LiFePO$_4$/graphite lithium-ion battery,'' \emph{J.
Power Sources}, vol.~195, no.~9, pp.~2961--2968, 2010,
doi: 10.1016/j.jpowsour.2009.10.105.""",
}

for path in (r"paper\main.tex", r"paper2\main.tex"):
    T = io.open(path, encoding="utf-8").read()
    n = 0
    for key, entry in VERIFIED.items():
        pat = re.compile(r"(\\bibitem\{" + key + r"\}\n)(.*?)(?=\n\n\\bibitem|\n\n\\end\{thebibliography\})",
                         re.S)
        if pat.search(T):
            T = pat.sub(lambda m: m.group(1) + entry, T)
            n += 1
    io.open(path, "w", encoding="utf-8").write(T)
    keys = re.findall(r"\\bibitem\{([^}]*)\}", T)
    print(f"{path}: {n} entries replaced, {len(keys)} bibitems -> {keys}")
