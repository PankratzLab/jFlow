package org.genvisis.common.gui;

import java.awt.GraphicsEnvironment;
import java.util.Objects;
import java.util.Properties;

import javax.activation.DataHandler;
import javax.activation.DataSource;
import javax.activation.FileDataSource;
import javax.mail.BodyPart;
import javax.mail.Message;
import javax.mail.Multipart;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeBodyPart;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.MimeMultipart;
import javax.swing.JOptionPane;

import org.pankratzlab.common.Logger;

/** This helper class is a catch-all for any exceptions that would not otherwise be caught. */
public final class ExceptionHandler implements Thread.UncaughtExceptionHandler {

  public static final String X11_ERROR_MSG_FORE = "Error occurred with X11 forwarding - ";
  public static final String X11_ERROR_DISABLED = "it's likely that X11 forwarding is disabled; please check your SSH client settings and try again.";
  public static final String X11_ERROR_XMING_REC = "it's likely that X11 forwarding is enabled but you are missing an X11 forwarding server (we recommend Xming - http://sourceforge.net/projects/xming/)";

  public Logger log;

  public void setLog(Logger log) {
    this.log = log;
  }

  @Override
  public void uncaughtException(Thread t, Throwable e) {
    // Notify user
    if (log != null) {
      log.reportError("Uncaught Exception in Thread {" + t.getName() + "}: " + e.getMessage());
      log.reportException(e);
    } else {
      System.err.println("Error - Uncaught Exception in Thread {" + t.getName() + "}:");
      e.printStackTrace();
    }
    // Report to devs via e-mail
    if (!GraphicsEnvironment.isHeadless()) {
      int result = JOptionPane.showConfirmDialog(null,
                                                 "Unexpected error occurred:\n" + e.getMessage()
                                                       + "\nWould you like to report this to the developers?",
                                                 "Unexpected Error", JOptionPane.YES_NO_OPTION,
                                                 JOptionPane.QUESTION_MESSAGE);
      if (result == JOptionPane.YES_OPTION) {
        // Recipient's email
        String genvisisErrors = "TODO";

        // Requires a gmail username and password to authenticate with gmail smtp
        // To make an app password: myaccount.gmail.com > Security tab > App Passwords
        String user = "TODO";
        String appPassword = "TODO";

        // Set required smtp properties
        // Note: port has changed over time
        // See https://developers.google.com/gmail/imap/imap-smtp
        Properties emailProperties = System.getProperties();
        emailProperties.put("mail.smtp.port", "587");
        emailProperties.put("mail.smtp.auth", "true");
        emailProperties.put("mail.smtp.starttls.enable", "true");

        // Get the default Session object.
        Session mailSession = Session.getDefaultInstance(emailProperties, null);
        try {
          // Create a default MimeMessage object.
          MimeMessage message = new MimeMessage(mailSession);

          // Set From: field
          message.setFrom(new InternetAddress(user));

          // Set To: field
          message.addRecipient(Message.RecipientType.TO, new InternetAddress(genvisisErrors));

          // Set Subject: field
          message.setSubject("Genvisis - Uncaught Exception Report");

          // Createa the actual message - just a stack trace
          String bodyText = e.getMessage();
          for (StackTraceElement ste : e.getStackTrace()) {
            bodyText += ("<br/>" + ste.toString());
          }

          // Create a multipart message
          Multipart multipart = new MimeMultipart();

          // Create the message part
          BodyPart messageBodyPart = new MimeBodyPart();

          // Now set the actual message
          messageBodyPart.setText(bodyText);

          // Set text message part
          multipart.addBodyPart(messageBodyPart);

          if (Objects.nonNull(log) && Objects.nonNull(log.getFilename())) {
            // Part two - attach log file
            messageBodyPart = new MimeBodyPart();
            String filename = log.getFilename();
            DataSource source = new FileDataSource(filename);
            messageBodyPart.setDataHandler(new DataHandler(source));
            messageBodyPart.setFileName(filename);
            multipart.addBodyPart(messageBodyPart);
          }

          message.setContent(multipart);

          // Send message
          String emailHost = "smtp.gmail.com";
          Transport transport = mailSession.getTransport("smtp");
          transport.connect(emailHost, user, appPassword);
          transport.sendMessage(message, message.getAllRecipients());
          transport.close();
        } catch (Exception exc) {
          exc.printStackTrace();
        }
      }
    }
  }
}
