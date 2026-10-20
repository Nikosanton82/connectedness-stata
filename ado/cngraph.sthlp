{smcl}
{* *! version 0.7.2  20oct2026}{...}
{title:Title}

{phang}
{bf:cngraph} {hline 2} Graph time-varying connectedness measures after connectedness


{title:Syntax}

{p 8 16 2}
{cmd:cngraph} {it:measure} [{cmd:,} {it:options}]

{p 4 4 2}
{it:measure} is one of

{synoptset 10}{...}
{synopt:{cmd:tci}}total connectedness index{p_end}
{synopt:{cmd:to}}total directional connectedness TO others, one panel per variable{p_end}
{synopt:{cmd:from}}total directional connectedness FROM others{p_end}
{synopt:{cmd:net}}net total directional connectedness{p_end}
{synopt:{cmd:npt}}net pairwise transmission count{p_end}
{synopt:{cmd:npdc}}net pairwise directional connectedness, one panel per pair{p_end}
{synopt:{cmd:pci}}pairwise connectedness index{p_end}
{synopt:{cmd:inf}}pairwise influence index{p_end}
{synopt:{cmd:gtci}}group TCIs, after {cmd:estat group, type(internal)} or {cmd:type(external)}{p_end}
{synopt:{cmd:tciext}}TCI_ext as defined in ConnectednessApproach, after {cmd:estat group, type(aggregated)}; see {helpb connectedness postestimation}{p_end}

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt fr:ame(name)}}graph the measures in frame {it:name}, for example the frame written by {cmd:estat group}; default is the results frame of the last estimates{p_end}
{synopt:{opt comp:are(namelist)}}overlay the same measures from stored {cmd:connectedness} results{p_end}
{synopt:{opt raw}}plot uncorrected TCIs even if {cmd:corrected} was specified{p_end}
{synopt:{opt sel:ect(name)}}show only panels involving {it:name}{p_end}
{synopt:{opt ban:d(#)}}graph frequency band {it:#} after {cmd:approach(frequency)}; default is the total{p_end}
{synopt:{opt comp:onent(name)}}graph the {cmd:contemporaneous} or {cmd:lagged} component after {cmd:approach(r2)}; default is the overall measure{p_end}
{synopt:{opt ti:tle(string)}}overall title{p_end}
{synopt:{opt name(name[, replace])}}name of the graph{p_end}
{synopt:{opt sav:ing(filename[, replace])}}save the graph{p_end}
{synopt:{opt nodraw}}do not display the graph{p_end}
{synopt:{opt combopts(string)}}options passed to {helpb graph combine}{p_end}
{synopt:{opt nosty:le}}use the colors and look of the current graph scheme{p_end}
{synopt:{it:twoway_options}}options passed to each panel{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:cngraph} draws the measures stored in a results frame of {helpb connectedness} or
{cmd:estat group}. It reads the variable names and settings from the frame's
characteristics.

{pstd}
Panel titles follow the R package: in {cmd:npdc}, a panel titled "A - B" shows the net
transmission from A to B, so positive values mean A is a net transmitter to B. More
than 36 panels require {opt select()}.

{pstd}
Overlays: with {opt compare()}, the measures of each stored estimate are drawn as
dashed lines, in other colors, over the current ones, matched on the time variable, so rolling-window and
TVP-VAR results can be compared; the legend labels the current results by model
(for example, TVP-VAR) and the overlays by the names of the stored estimates. After {cmd:connectedness, asymmetric}, the measures of
positive and negative changes are overlaid automatically. Overlaid panels are drawn as
lines and a note gives the line patterns; single panels show areas, as in the R
package.

{pstd}
After {cmd:connectedness, approach(frequency)}, {cmd:cngraph tci} draws the total
TCI and the TCI of every band in one graph; the other measures show the total
unless {opt band(#)} picks a band. After {cmd:approach(r2)} with lags,
{cmd:cngraph tci} draws the overall, contemporaneous and lagged TCIs, and
{opt component()} picks a component for the other measures. After {cmd:approach(joint)}, only {cmd:tci},
{cmd:to}, {cmd:from} and {cmd:net} are available.

{pstd}
The graphs share the look of {helpb cnnetwork}: a light background, muted axes and
gridlines, and a fixed palette (dark blue, mustard, teal, rust, grey and plum) for the
series of one graph. NET and NPDC panels are drawn as areas, dark blue above zero (net
transmitter) and mustard below zero (net receiver), the colors of net transmitters and
receivers in {helpb cnnetwork}; TO and FROM panels are dark blue areas. Overlaid results
are drawn in the next colors of the palette and with dashed patterns. {opt nostyle}
keeps the colors and look of the current graph scheme instead; {it:twoway_options} are
applied last and can change any element.

{pstd}
The horizontal axis shows the dates of the results frame (see the {opt datevar()}
option of {helpb connectedness}); daily, monthly and quarterly dates are labeled by
year.


{title:Examples}

{phang2}{cmd:. use cn_acg2020, clear}{p_end}
{phang2}{cmd:. tsset month}{p_end}
{phang2}{cmd:. connectedness EUR GBP JPY CHF, horizon(12) window(200) frame(roll) replace}{p_end}
{phang2}{cmd:. estimates store roll200}{p_end}
{phang2}{cmd:. connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) frame(tvp) replace}{p_end}
{phang2}{cmd:. cngraph tci, compare(roll200)}{p_end}
{phang2}{cmd:. cngraph net}{p_end}
{phang2}{cmd:. cngraph npdc, select(EUR)}{p_end}

{phang2}{cmd:. use cn_gg2018, clear}{p_end}
{phang2}{cmd:. tsset month}{p_end}
{phang2}{cmd:. connectedness US_FPU-JP_CPU, model(tvpvar) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. estat group, type(internal) groups("US: US_FPU US_MPU US_TPU US_CPU" "JP: JP_FPU JP_MPU JP_TPU JP_CPU") replace}{p_end}
{phang2}{cmd:. cngraph gtci, frame(cn_group)}{p_end}
{phang2}{cmd:. cngraph net, frame(cn_group)}{p_end}


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
Help:  {helpb connectedness}, {helpb connectedness postestimation}, {helpb cnnetwork}{p_end}
