{smcl}
{* *! version 0.2.1  18oct2026}{...}
{title:Title}

{phang}
{bf:cnnetwork} {hline 2} Network plot of pairwise connectedness after connectedness


{title:Syntax}

{p 8 16 2}
{cmd:cnnetwork} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt m:ethod(npdc|pci)}}net pairwise directional connectedness (default) or
pairwise connectedness index{p_end}
{synopt:{opt thr:eshold(#)}}omit edges whose normalized strength is at most {it:#}, 0 <= {it:#} < 1; default 0{p_end}
{synopt:{opt at(exp)}}draw one period, matched on the time or date variable; default is the average over all periods{p_end}
{synopt:{opt fr:ame(name)}}use the measures in frame {it:name}; default is the results frame of the last estimates{p_end}
{synopt:{opt ti:tle(string)}}title{p_end}
{synopt:{opt name(name[, replace])}}name of the graph{p_end}
{synopt:{opt sav:ing(filename[, replace])}}save the graph{p_end}
{synopt:{opt nodraw}}do not display the graph{p_end}
{synopt:{it:twoway_options}}other options passed to {helpb twoway}{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:cnnetwork} draws the network of pairwise connectedness, following the R function
{cmd:PlotNetwork()} (Gabauer 2022). Variables sit on a circle. With {cmd:method(npdc)},
an arrow runs from the net transmitter to the net receiver of each pair, its width
proportional to the net pairwise directional connectedness normalized by the largest
value; net transmitters are drawn as filled circles and net receivers as hollow circles,
sized by the absolute value of NET. Colors follow the current graph scheme.
With {cmd:method(pci)}, undirected edges have widths proportional to the pairwise
connectedness index. Edge widths are drawn in five classes.

{pstd}
{cmd:method(npdc)} also works on frames written by {cmd:estat group}, for example to
draw the network of aggregated groups.


{title:Examples}

{phang2}{cmd:. use cn_dy2012, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. cnnetwork}{p_end}
{phang2}{cmd:. cnnetwork, method(pci) threshold(0.2)}{p_end}
{phang2}{cmd:. cnnetwork, at(td(15sep2008))}{p_end}


{marker author}{...}
{title:Author}

{pstd}
Nikolaos Antonakakis{break}
Department of Accounting, Economics and Finance{break}
University of Nicosia, UNIC Athens{break}
Athens, Greece{break}
{browse "mailto:antonakakis.n@unic.ac.cy":antonakakis.n@unic.ac.cy}

{pstd}
The numerical routines port the {it:ConnectednessApproach} R package by David
Gabauer, distributed under the GNU General Public License version 3; this package
is distributed under the same license.


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume ??, number ??: {browse "https://doi.org/":st0001}{p_end}

{p 7 14 2}
Help:  {helpb connectedness}, {helpb connectedness postestimation}, {helpb cngraph}{p_end}
