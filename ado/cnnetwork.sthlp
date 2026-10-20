{smcl}
{* *! version 0.7.1  20oct2026}{...}
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
{synopt:{opt ban:d(#)}}draw frequency band {it:#} after {cmd:approach(frequency)}; default is the total{p_end}
{synopt:{opt comp:onent(name)}}draw the {cmd:contemporaneous} or {cmd:lagged} component after {cmd:approach(r2)}{p_end}
{synopt:{opt str:aight}}straight edges instead of curved ones{p_end}
{synopt:{opt curv:ature(#)}}bend of the edges, from -1 to 1; default 0.2{p_end}
{synopt:{opt col:ors(color1 color2)}}colors of net transmitters and net receivers{p_end}
{synopt:{opt lab:els(where)}}{cmd:auto} (default), {cmd:inside} or {cmd:outside} the nodes{p_end}
{synopt:{opt nonote}}omit the note on edge widths{p_end}
{synopt:{opt ti:tle(string)}}title{p_end}
{synopt:{opt name(name[, replace])}}name of the graph{p_end}
{synopt:{opt sav:ing(filename[, replace])}}save the graph{p_end}
{synopt:{opt nodraw}}do not display the graph{p_end}
{synopt:{it:twoway_options}}other options passed to {helpb twoway}{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:cnnetwork} draws the network of pairwise connectedness, following the R function
{cmd:PlotNetwork()} (Gabauer 2022). Variables sit on a circle, starting at the top and
running clockwise in {it:varlist} order. With {cmd:method(npdc)}, a curved arrow runs
from the net transmitter to the net receiver of each pair; its width and shade grow
with the net pairwise directional connectedness, normalized by the largest value and
drawn in six classes. Net transmitters are dark blue and net receivers mustard, and
nodes are slightly larger for larger absolute NET; a legend below the circle explains
the colors. With {cmd:method(pci)}, undirected edges have widths and shades that grow
with the pairwise connectedness index. The subtitle gives the period: the average
over the first to the last date of the results frame, or the date in {opt at()}.

{pstd}
The graph has its own light background and does not depend on the graph scheme.
Nodes, edges and labels shrink as the number of variables grows. With up to 12
variables whose names have at most five characters, the names are written inside the
nodes; otherwise they are placed outside the circle. {it:twoway_options} are applied
last and can change, for example, the title or the graph size.

{pstd}
{cmd:method(npdc)} also works on frames written by {cmd:estat group}, for example to
draw the network of aggregated groups.

{pstd}
After {cmd:connectedness, approach(frequency)}, {opt band(#)} draws the network of
one frequency band. {cmd:cnnetwork} is not available after {cmd:approach(joint)},
which has no pairwise measures; {cmd:approach(extjoint)} has them.


{title:Options}

{phang}
{opt straight} draws straight edges. {opt curvature(#)} sets how far the middle of
each edge bends away from the straight line, as a share of its length; negative
values bend the other way.

{phang}
{opt colors(color1 color2)} sets the colors of net transmitters and net receivers,
for example {cmd:colors(navy maroon)} or {cmd:colors("44 95 125" "201 154 46")} (the
default). With {cmd:method(pci)} all nodes take {it:color1}.

{phang}
{opt labels(inside|outside)} overrides the automatic placement of the variable names;
long names written inside small nodes may overflow them.

{phang}
{opt nonote} omits the note on edge widths and the threshold.


{title:Examples}

{phang2}{cmd:. use cn_dy2012, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. cnnetwork}{p_end}
{phang2}{cmd:. cnnetwork, method(pci) threshold(0.2)}{p_end}
{phang2}{cmd:. cnnetwork, at(td(15sep2008))}{p_end}
{phang2}{cmd:. cnnetwork, straight colors(navy maroon)}{p_end}


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
