<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>prova</title>
 </head>
 <body >

<?PHP

if (!get_magic_quotes_gpc()) {
	$email_iscritto = addslashes($_REQUEST["email_iscritto"]);
} else {
	$email_iscritto = $_REQUEST["email_iscritto"];
}


?>

<BR>

<?PHP

require("./class.phpmailer.php");

$mail = new PHPMailer();
$mail->IsSMTP();
$mail->Host = "localhost";
$mail->From = "tdstool@unimore.it"; // indirizzo email del mittente (possibilmente da cambiare)
$mail->FromName = " NUOVO ISCRITTO "; // nome e cognome del mittente


// metto nella variabile txtSubject l'oggetto del messaggio
$txtSubject = "Oggetto del messaggio";

// metto nella variabile txtBody il testo del messaggio in formato html (si puo fare anche un email solo in testo semplice ascii)

$txtBody = "
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" 
 \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">
<html><body>
Buon giorno.<br>
<b>Questo &egrave; un messaggio di prova.</b><br>
email del mittente: $email_iscritto<br>
</body></html>
";

// metto nella variabile txtAltBody il messaggio in testo semplice ascii
$txtAltBody  = "
buon giorno\n
questo e' un messaggio di prova\n
email del mittente: $email_iscritto\n
";
   
$mail->AddAddress("tdstool@unimore.it", "Davide Calanca");  // da cambiare, mi raccomando ... 
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


