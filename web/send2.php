<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>prova</title>
 </head>
 <body >

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
$mail->FromName = " NUOVO ISCRITTO "; // nome e cognome del mittente


// metto nella variabile txtSubject l'oggetto del messaggio
$txtSubject = "Newsletter Subscription";

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
   
$mail->AddAddress("tdstool@s3.infm.it", "TDS Tool Subscription");  // da cambiare, mi raccomando ... 
$mail->IsHTML(true);   // mettere false se il messaggio e' in testo semplice
$mail->SMTPAuth = false;
//$mail->Username = "smtp.fisica.web";
//$mail->Password = "ruphecrydy";
$mail->Subject = $txtSubject;
$mail->Body    = $txtBody;
$mail->AltBody = $txtAltBody;

if(!$mail->Send())
{
	die("errore: ".$mail->ErrorInfo);
} else {
	print("iscrizione avvenuta con successo<br>");
}
   
$mail->ClearAllRecipients();
?>


</body>
</html>


