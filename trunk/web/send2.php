<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
<head>
 <title>TDStool Newsletter Subscription</title>
 <link rel="stylesheet" href="style.css" type="text/css" media="screen" />
</head>

<body >
<div class="PageBackgroundSimpleGradient"> </div>

<div class="Main">
 <div class="Sheet">
  <div class="Sheet-tl"></div>
  <div class="Sheet-tr"><div></div></div>
  <div class="Sheet-bl"><div></div></div>
  <div class="Sheet-br"><div></div></div>
  <div class="Sheet-tc"><div></div></div>
  <div class="Sheet-bc"><div></div></div>
  <div class="Sheet-cl"><div></div></div>
  <div class="Sheet-cr"><div></div></div>
  <div class="Sheet-cc"></div>
  <div class="Sheet-body">
   <div class="Header">
    <div class="Header-png"></div>
    <div class="Header-jpeg"></div>
    <div class="logo"></div>
   </div>
   
   <div class="contentLayout">
    <div class="content">
     <div class="Post">
      <div class="Post-body">
       <div class="Post-inner">
   

<?PHP

if (!get_magic_quotes_gpc()) {
	$email_iscritto = addslashes($_REQUEST["email_iscritto"]);
	$nome_iscritto = addslashes($_REQUEST["nome_iscritto"]);
	$cognome_iscritto = addslashes($_REQUEST["cognome_iscritto"]);
	$affiliazione_iscritto = addslashes($_REQUEST["affiliazione_iscritto"]);
} else {
	$email_iscritto = $_REQUEST["email_iscritto"];
	$nome_iscritto = $_REQUEST["nome_iscritto"];
	$cognome_iscritto = $_REQUEST["cognome_iscritto"];
	$affiliazione_iscritto = $_REQUEST["affiliazione_iscritto"];
}


?>

<BR>

<?PHP

require("./class.phpmailer.php");

$mail = new PHPMailer();
$mail->IsSMTP();
$mail->Host = "localhost";
$mail->From = "tdstool@s3.infm.it"; // indirizzo email del mittente (possibilmente da cambiare)
$mail->FromName = "Tdstool website"; // nome e cognome del mittente


// metto nella variabile txtSubject l'oggetto del messaggio
$txtSubject = "Tdstool Newsletter NUOVO ISCRITTO";

// metto nella variabile txtBody il testo del messaggio in formato html (si puo fare anche un email solo in testo semplice ascii)

$txtBody = "
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" 
 \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">
<html><body>
email del mittente: $email_iscritto<br>
Nome: $nome_iscritto<br>
Cognome: $cognome_iscritto<br>
Affiliazione: $affiliazione_iscritto<br>
</body></html>
";

// metto nella variabile txtAltBody il messaggio in testo semplice ascii
$txtAltBody  = "
email del mittente: $email_iscritto\n
Nome: $nome_iscritto\n
Cognome: $cognome_iscritto\n
Affiliazione: $affiliazione_iscritto\n
";
   
$mail->AddAddress("tdstool@s3.infm.it", "TDS Tool Subscription");
$mail->IsHTML(true);   // mettere false se il messaggio e' in testo semplice
$mail->SMTPAuth = false;
//$mail->Username = "smtp.fisica.web";
//$mail->Password = "ruphecrydy";
$mail->Subject = $txtSubject;
$mail->Body    = $txtBody;
$mail->AltBody = $txtAltBody;

if(!$mail->Send())
{
	die("SUBSCRIPTION ERROR (Please contact tdstool@s3.infm.it). Error code: ".$mail->ErrorInfo);
} else {
	print("<span>");
	print("<b>Newsletter subscription successfull !</b></br></br>");
	print("</span>");
}
   
$mail->ClearAllRecipients();
?>

<a href="index.php"><span><span>Back to Home</span></span></a>


       </div>
      </div>
     </div>
    </div>
   </div>
  </div>
 </div>
</div>

</body>
</html>


