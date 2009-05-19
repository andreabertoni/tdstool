<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" dir="ltr" lang="en-US"><head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<meta http-equiv="X-UA-Compatible" content="IE=EmulateIE7" /><title>TDStool Subscribe</title>

<script type="text/javascript" src="script.js"></script>
<link rel="stylesheet" href="style.css" type="text/css" media="screen" />
<!--[if IE 6]><link rel="stylesheet" href="style.ie6.css" type="text/css" media="screen" /><![endif]-->
</head>
<body>
<div class="PageBackgroundSimpleGradient"> </div>
<div class="Main">
<div class="Sheet">
<div class="Sheet-tl"></div>
<div class="Sheet-tr">
<div></div>
</div>
<div class="Sheet-bl">
<div></div>
</div>
<div class="Sheet-br">
<div></div>
</div>
<div class="Sheet-tc">
<div></div>
</div>
<div class="Sheet-bc">
<div></div>
</div>
<div class="Sheet-cl">
<div></div>
</div>
<div class="Sheet-cr">
<div></div>
</div>
<div class="Sheet-cc"></div>
<div class="Sheet-body">
<div class="Header">
<div class="Header-png"></div>
<div class="Header-jpeg"></div>
<div class="logo">
</div>
</div>

<?php include("navbar_up.php"); ?>

<div class="contentLayout">
<div class="content">
<div class="Post">
<div class="Post-body">
<div class="Post-inner">
<h2 class="PostHeaderIcon-wrapper"><img src="images/PostHeaderIcon.png" alt="PostHeaderIcon" height="12" width="12" />
<span class="PostHeader">Subscribe to TDStool Newsletter</span>
</h2>
<br />

<form method="post" enctype="multipart/form-data" action="send2.php">
<table width="100%">
  <tr>
    <td width="27%">E-mail:<br />
        <span class="txtlittle"></span> <br />
      </td>
    <td  width="73%"><input type="text" id="Text1" name="email_iscritto" size="36" /> <br /><br />
    </td>
  </tr>
  <tr>
    <td width="27%">Optional fields: <br />
        <span class="txtlittle"></span> <br />
      </td>
    <td  width="73%"> <br />
    </td>
  </tr>
  <tr>
    <td width="27%">First Name:</td>
    <td  width="73%"><input type="text" name="nome_iscritto" id="FirstName4" size="36" /></td>
  </tr>
  <tr>
    <td width="27%">Last Name:</td>
    <td  width="73%"><input type="text" name="cognome_iscritto" id="LastName2" size="36" /></td>
  </tr>
  <tr>
    <td width="27%">Affiliation:<br />
      <span class="txtlittle">(Department; Institute; University; other)</span> <br /> </td>
    <td  width="73%">
      <input type="text" name="affiliazione_iscritto" id="affiliation" size="36" /> </td>
  </tr>
  <tr>
    <td colspan="2" class="centerTab"><button class="Button" type="submit" value="iscrizione">
<span class="btn"><span class="t">Subscribe</span> <span class="r"><span></span></span>
<span class="l"></span> </span> </button></td>
  </tr>
</table>
</form>

</div>
</div>
</div>
</div>
<div class="sidebar1">

<?php include("highlights_box.php"); ?>
<!--?php include("newsletter_box.php"); ?-->
<?php include("contact_box.php"); ?>

</div>
</div>
</div>
<div class="cleared"></div>

<?php include("footer_dn.php"); ?>

</div>
</div>
<div class="cleared"></div>
</body></html>