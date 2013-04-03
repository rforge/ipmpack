
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
            <table width="1045px" style="table-layout:fixed;background:#FFFFFF" cellpadding="10">
              <tbody>
                <tr>
                  <td valign="top" width="350">
                    Under construction<p>

For performing bayesian IPM analyses the following function might be of use:<p>  
## A function that outer can use showing numbers from x to y via production, growth, survival and distribution offspring<p>
.fecPostCensus <- function(x,y,cov=data.frame(covariate=1),fecObj, growObj,survObj) {<br>
	newd <- data.frame(cbind(cov,size=x),
			stringsAsFactors = FALSE)<br>
	newd$size2 <- x^2<br>
	newd$size3 <- x^3<br>
	if (length(grep("expsize",fecObj@offspringRel$formula))>0 |
			length(grep("expsize", growObj@fit$formula))>0) { newd$expsize <- exp(x)}<br>            
	if (length(grep("logsize",fecObj@offspringRel$formula))>0 |
			length(grep("logsize", growObj@fit$formula))>0) { newd$logsize <- log(x)}<br>            
	u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
			dnorm(y,predict(fecObj@offspringRel, newdata=newd, type="response"), fecObj@sdOffspringSize) *
			surv(size=x, cov=cov, survObj=survObj)<br>
	return(u)<br>
}<br>

                    
                  </td>
                </tr>
              </tbody>
            </table>
    <br><br><br><br><br><br><br><br><br><br><br>
  </body>
</html>
