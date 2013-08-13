
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';


?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>IPMpack</title>
    <link rel="stylesheet" type="text/css" href="ipmpackstyle.css" /> 
	<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-4964902-4']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

  </script>
  </head>

  <body style="background:#CCCCCC">
    <a name="top"></a>
    <table style="vertical-align:middle;background:#CCCCCC" width="100%" height="100%">
      <tbody>
        <tr align="center">
          <td>
            <table style="vertical-align:middle;background:#1F1209" width="1055px">
              <tbody><tr align="center"><td>
            <!-- 1.- HEADER:-->
            <table style="background-color:#47697E" width="1050"  style="table-layout:fixed">
              <tbody>
                <tr>
                  <td width="10px">&nbsp;</td>
                  <td width="20px">&nbsp;</td>
                  <td width="920px" style="vertical-align:middle">
                    <p style="white-space:nowrap;vertical-align:middle;margin-top:0px;margin-bottom:0px;font-size:38px;letter-spacing:0.25em;color:#FFCC33">IPMpack</p>
                    <p style="white-space:nowrap;vertical-align:middle;font-size:24px;color:#FCF1D1">an R package for Integral Projection Models</p>
                  </td>
                  <td valign="middle" width="100px"><img alt="IPMpack logo" border="0" height="100px" src="IPMpacklogo.jpg?width=100px"></td>
                </tr>
              </tbody>
            </table>
            <!-- 2.- LINKS:-->
            <table width="1050px" style="table-layout:fixed;vertical-align:middle" cellspacing="5px">
              <tbody>
                <tr>
                  <td class="button"><a href="#summary" class="button">SUMMARY</a></td>
                  <td class="button"><a href="#description" class="button">PACKAGE</a></td>
                  <td class="button"><a href="#refs" class="button">REFERENCES</a></td>
                  <td class="button"><a href="#faqs" class="button">FAQ</a></td>
                </tr>
              </tbody>
            </table>
            <!-- 3.- LOCAL LINKS:-->
            <table width="1045px" style="table-layout:fixed;background:#FFFFFF" cellpadding="10">
              <tbody>
                <tr>
                  <td align="center" valign="top" width="350">
                    <p><br></p>
                    <img alt="IPM'd species" border="0" width="350px" src="IPMpackPhotos.jpg">
                    <p align="left"><br>&nbsp<a href="http://r-forge.r-project.org/"><img height="30px" src="Rforgelogo.png" border="0" alt="R-Forge Logo"></a></p>
                  </td>
                  <td width="400" valign="top">
                    <h1><b>IPMpack</b></h1>
                    <p style="font-size:16px;line-height:1.25"><b>Authors:</b><br><br><a href="http://www.zoo.ox.ac.uk/people/view/metcalf_cj.htm" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">C. Jessica E. Metcalf</a> University of Oxford<br><a href="http://www.stri.si.edu/english/scientific_staff/staff_scientist/scientist.php?id=68" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">Sean M. McMahon</a> Smithsonian Tropical Research Institute<br><a href="https://sites.google.com/site/robresearchsite/" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">Roberto Salguero-Gómez</a> MPIDR<br><a href="http://www.ru.nl/animalecology/staff/jongejans/" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">Eelke Jongejans</a> Radboud University Nijmegen<br><a href="http://www.eeb.uconn.edu/people/cmerow/home/Home.html" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">Cory Merow</a> STRI and University of Connecticut<br></p>
                    <p style="font-size:16px;line-height:1.25;text-align:justify"><br><b>Developed at:</b><br><br><a href="http://www.demogr.mpg.de" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">Max Planck Institute for Demographic Research</a></p>
                    <p style="font-size:16px;line-height:1.25;text-align:justify"><br><b>Cite as:</b><br><br>Metcalf CJE, McMahon SM, Salguero-Gómez R & Jongejans E (2013) IPMpack: an R package for Integral Projection Models. <i>Methods in Ecology and Evolution</i> 4:195-200 <a href="www.ru.nl/publish/pages/654854/metcalf_2013_methods_ecol_evol.pdf">PDF</a></p>
                    <p style="font-size:16px;line-height:1.25;text-align:justify"><br><b>Resources:</b><br><br>Subscribe to the <a href="http://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/ipmpack-users" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">users email-list</a><br><br>Download the R package <a href="http://cran.r-project.org/web/packages/IPMpack/index.html" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">IPMpack</a> (latest version: 2.0 (Jul'13). In v2.0 we have substantially streamlined IPMpack by homogenizing function naming and argument conventions and consolidating related functions in an effort to make functions easier to find and remember. Due to name changes, some of your past IPMpack code may not function, though we've made every effort to retain functions temporarily to ease the transition to version 2.0. Some deprecated functions are now hidden, which means that they can be accessed by preceeding the function name with a '.' (e.g. .generateDataDiscrete()). Don't hesitate to contact us through the <a href="http://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/ipmpack-users" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">list serv</a> to ask for tips on updating your code to 2.0. Our new naming system is now stable, unless catastrophe strikes. Please update the function names in your R code as indicated in the warning messages, because these functions will be discontinued in future versions of IPMpack.</p>
                  </td>
                  <td width="50px"></td>
                </tr>
              </tbody>
            </table>
            <!-- 4.- SECTIONS:-->
            <!-- 4.1- Summary:-->
            <a name="summary"></a>
            <table  valign="middle" width="1050px" style="table-layout:fixed;background:#1F1209">
              <tbody>
                <tr class="title"><td colspan="2" class="section">Project summary</td></tr>
                <tr style="background:#FFFFFF">
                  <td class="main">
                    <p class="parag"><b>IPMpack</b> is an R package (R Development Core Team 2013) that allows users to build and analyse Integral Projection Models. An IPM is a demographic tool to explore the dynamics of populations where individuals' fates depend on state variables that are continuous (e.g. weight, diameter at breast height, height, limb length, rosette diameter...) or quasi-continuous (number of leaves, age, number of reproductive structures) and may be a mixture of discrete (e.g. seedbank) and continuous. IPMs track the distribution of individuals <i>n</i> across these state variables between census times (e.g., year <i>t</i> and year <i>t+1</i>) by projecting from models that define the underlying vital rates (e.g., survival, growth, reproduction) as a function of the (quasi-)continuous state variables. Version 1.5 of <a href="http://cran.r-project.org/web/packages/IPMpack/index.html" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">IPMpack</a> is now available on CRAN /  R-Forge. For those who wish to try it, it can be installed by opening the R console and typing: <br></p>
										<p><code>install.packages("IPMpack", repos = "http://R-Forge.R-project.org", type = "source")</code>
					<p><br></p>

                    <a href="#top" class="totop">Back to top</a>
                  </td>
                  <td class="main">
                    <p><a href="file:///C:/Users/Eelke/Documents/Eclipse/IPMpack/www/fig1.jpg?width=400px" target="blank"><img style="margin:0px;padding:0px;border:none" width="420px" align="middle" src="fig1.jpg?width=400px"></a><br><br></p>
                    <p class="caption">Fig. 1. <b>IPMpack</b> output examples.</p>
										<p class="parag">We have set up a IPMpack <a href="http://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/ipmpack-users" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">users email-list</a> so users can ask questions or provide comments, suggestions or criticism that can help us improve the package. <br></p>
                  </td>
                </tr>
              </tbody>
            </table>
            <!-- 4.2- Package:-->
            <a name="description"></a>
            <table  valign="middle" width="1050px" style="table-layout:fixed;background:#1F1209">
              <tbody>
                <tr class="title"><td colspan="2" class="section">Package description</td></tr>
                <tr style="background:#FFFFFF">
                  <td class="main">
                    <p class="parag">To use <b>IPMpack</b>'s full capacities, it is helpful if the data are in a specific format in R, i.e., a dataframe with the following variables and column names, where each
row represents one measurement in the population:</p>
                    <ul>
                      <li><b>size</b> of individuals in census time <i>t</i>.</li>
                      <li><b>sizeNext</b> of individuals in census time <i>t+1</i>.</li>
                      <li><b>surv</b>: survival of individuals from census time <i>t</i> to <i>t+1</i> (a 0 or 1).</li>
                      <li><b>fec1,...</b>: as many columns as desired relating size to sexual reproduction.</li>
                      <li><b>stage</b> of individuals in census time <i>t</i> needs to be specified if you want to include discrete classes. For rows in the dataframe where <b>size</b> is not an NA, this must be the word <i>continuous</i>.</li>
                      <li><b>stageNext</b> of individuals in census time <i>t+1</i> (similar to <b>stage</b>).</li>
                      <li><b>number</b> of individuals corresponding to each row in the dataframe.</li>
                      <li><b>covariate</b>: value of a discrete covariate in census time <i>t</i>, such as light environment, age or patch type.</li>
                      <li><b>covariateNext</b>: value of a discrete covariate in census time <i>t+1</i>.</li>
                      <li>...any other covariates of interest, named as desired (precipition, habitat, temperature,...).</li>
                    </ul>
                    <p><br></p>
                    <a href="#top" class="totop">Back to top</a>
                  </td>
                  <td style="vertical-align:top;padding:40px" width="420px">
                    <a href="fig2.jpg" target="blank"><img style="margin:0px;padding:0px;border:none;width:420px" align="middle" src="fig2.jpg?width=420px"></a>
                    <p class="caption">Fig. 2. <b>IPMpack</b> usage flow-chart.</p><br>
                    <ul>
                      <li>Fits a range of vital rate functions.</li>
                      <li>Constructs survival-growth and fecundity IPMs.</li>
                      <li>Provides diagnostics for contructed IPMs (see Fig1 above).</li>
                      <li>Offers a large and growing number of IPM analysis tools.</li>
                    </ul>
                    <p style="text-align:justify;line-height:1.75;font-size:18px">Future versions will include:</p>
                    <ul style="text-align:justify;line-height:1.75;font-size:18px">
                      <li>Clonality IPMs.</li>
                      <li>Multiple continuous state variables in the same IPM.</li>
                      <li>etc. (suggestions are very welcome)</li>
                    </ul>
                  </td>
                </tr>
              </tbody>
            </table>
            <!-- 4.4- References:-->
            <a name="refs"></a>
            <table  valign="middle" width="1050px" style="table-layout:fixed;background:#1F1209">
              <tbody>
                <tr class="title"><td colspan="2" class="section">References</td></tr>
                <tr style="background:#FFFFFF">
                  <td style="vertical-align:middle;padding:40px" width="420px">
				    Papers consulted when we built IPMpack:
                    <p class="Refs">Caswell  H (2001) <b>Matrix Population Models: Analysis, Construction and Interpretation</b>. 2nd ed, Sinauer, Sunderland, Massachusetts<br><br></p>
					<p class="Refs">Childs DZ, Rees M, Rose KE, Grubb PJ & Ellner SP (2004) <b>Evolution of size-dependent flowering in a variable environment: Construction and analysis of a stochastic integral projection model.</b> <i>Proc Roy Soc B</i> 271:471–475 <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1691612/pdf/15101702.pdf" rel="nofollow" style="color:#84002E" target="_blank">pdf</a><br><br></p>
					<p class="Refs">Cochran ME & Ellner SP (1992) <b>Simple methods for calculating age-based life history parameters for stage-structured populations.</b> <i>Ecol Monogr</i> 62:345–364 <a href="http://dx.doi.org/10.2307/2937115" rel="nofollow" style="color:#84002E" target="_blank">doi</a><br><br></p>
					<p class="Refs">Ellner SP & Rees M (2006) <b>Integral projection models for species with complex life-histories.</b> <i>Am Nat</i> 167:410-428<br><br></p>
					<p class="Refs">Metcalf CJE, Horvitz CC, Tuljapurkar S & Clark DA (2009) <b>A time to grow and a time to die: a  new way to analyze the dynamics of size, light, age and death of tropical trees.</b> <i>Ecology</i> 90:2766–2778 <a href="http://dx.doi.org/10.1890/08-1645.1" rel="nofollow" style="color:#84002E" target="_blank">doi</a><br><br></p>
					<p class="Refs">Rees M & Rose KE (2002) <b>Evolution of flowering strategies in <i>Oenothera glazioviana</i>: an integral projection model approach.</b> <i>Proc Roy Soc B</i> 269:1509-1515 <a href="http://rspb.royalsocietypublishing.org/content/269/1499/1509.full.pdf" rel="nofollow" style="color:#84002E" target="_blank">pdf</a><br><br></p>
                    <p class="Refs">Tuljapurkar S (1990) <b>Population Dynamics in Variable Environments</b>. Springer, Berlin<br><br></p>
					<p class="Refs">Zuidema PA, Jongejans E, Chien PD, During HJ & Schieving F (2010) <b>Integral Projection Models for trees: a new parameterization method and a validation of model output.</b> <i>Journal of Ecology</i> 98:345-355 <a href="www.ru.nl/publish/pages/545422/zuidema_et_al_-_2010_-_journal_of_ecology.pdf" rel="nofollow" style="color:#84002E" target="_blank">pdf</a>
                    <p><br></p>
                    <a href="#top" class="totop">Back to top</a>
                  </td>
                  <td style="vertical-align:top;padding:40px" width="420px">
				    Papers using IPMpack:
                    <p class="Refs">Metcalf CJE, McMahon SM, Salguero-Gómez R & Jongejans E (2013) <b>IPMpack: an R package for Integral Projection Models</b>. <i>Meth Ecol Evol</i> 4:195-200 <a href="www.ru.nl/publish/pages/654854/metcalf_2013_methods_ecol_evol.pdf">PDF</a><br><br></p>
                    <p class="Refs">Dauer JT & Jongejans E (2013) <b>Elucidating the population dynamics of Japanese knotweed using integral projection models</b>. <i>PLoS ONE</i> in press<br><br></p>
					<p><br><br><br></p>
					<p><br><br><br></p>
					<p>A list of all known IPM papers can be downloaded here: <a href="IPMpublications.xlsx">IPMpublications</a>. Please <a href="mailto:E.Jongejans@science.ru.nl?subject=newOrMissedIPMPaper">email us</a> if we missed (a new) one.
					<p><br><br><br></p>
					<p>Photo credits: <br><i>Xanthoparmelia</i> lichens by Ann Pringle<br>crocs by Owen Jones<br><i>Succisa pratensis</i> by Lidewij Keser
					<p>Webdesign credits: <br>
					Fernando Colchero of <a href="http://basta.r-forge.r-project.org/" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">BaSTA</a>
                  </td>
                </tr>
              </tbody>
            </table>
            <!-- 4.4- References:-->
            <a name="faqs"></a>
            <table  valign="middle" width="1050px" style="table-layout:fixed;background:#1F1209">
              <tbody>
                <tr class="title"><td colspan="2" class="section">FAQs and other pointers</td></tr>
                <tr style="background:#FFFFFF">
                  <td style="vertical-align:middle;padding:40px" width="420px">
				    Here we list frequently asked questions, along with other pointers that we deem useful. Some of the answers are accompanied with links to R code on separate pages. If you have additional questions, please subscribe to the <a href="http://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/ipmpack-users" rel="nofollow" style="color:#84002E;font-family:verdana" target="_blank">users email-list</a> and post your questions there. 
                    <p class="Refs">Can I perform <b>bayesian</b> analyses with <i>IPMpack</i>? <p>No, we do not include functions for bayesian analyses. However, some R code (from previous versions of <i>IPMpack</i>) can be found <a href="bayes.php">here</a>.
                    <p class="Refs">How do I <b>plot</b> my matrices? <p>From version 2.0 <i>IPMpack</i> no longer contains the function countourPlot. Instead we recommend you use image.plot of the 'fields' package.
                    <p class="Refs">How do I compare my IPMs with <b>matrix</b> models?<p>The relationships stored in your vital rate objects can of course be used to parameterize matrix models as well as IPMs. On a separate <a href="bayes.php">page</a> we provide some code to get you started building matrix models from the vital rate objects yourself (we decided to focus our efforts on good-running IPM code and not to continue developing code for matrix models within <i>IPMpack</i>.
                    <p><br></p>
                    <a href="#top" class="totop">Back to top</a>
                  </td>
                  <td style="vertical-align:top;padding:40px" width="420px">
                  </td>
                </tr>
              </tbody>
            </table>
            <p><br><br></p>
          </td>
        </tr>
      </tbody>
    </table>
    <br><br><br><br><br><br><br><br><br><br><br>
    <div id="clustrmaps-widget"></div><script type="text/javascript">var _clustrmaps = {'url' : 'http://ipmpack.r-forge.r-project.org/', 'user' : 1004170, 'server' : '3', 'id' : 'clustrmaps-widget', 'version' : 1, 'date' : '2012-04-10', 'lang' : 'en', 'corners' : 'square' };(function (){ var s = document.createElement('script'); s.type = 'text/javascript'; s.async = true; s.src = 'http://www3.clustrmaps.com/counter/map.js'; var x = document.getElementsByTagName('script')[0]; x.parentNode.insertBefore(s, x);})();</script><noscript><a href="http://www3.clustrmaps.com/user/f33f528a"><img src="http://www3.clustrmaps.com/stats/maps-no_clusters/ipmpack.r-forge.r-project.org--thumb.jpg" alt="Locations of visitors to this page" /></a></noscript>
  </body>
</html>
