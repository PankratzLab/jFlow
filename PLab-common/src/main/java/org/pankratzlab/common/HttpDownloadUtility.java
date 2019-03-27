package org.pankratzlab.common;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import javax.xml.ws.http.HTTPException;

public class HttpDownloadUtility {

  public static class HttpDownloadException extends Exception {

    private static final long serialVersionUID = 1L;

    public HttpDownloadException(String message) {
      super(message);
    }

    public HttpDownloadException(Exception cause) {
      super(cause);
    }
  }

  private static final int BUFFER_SIZE = 4096;

  public static boolean canDownload(String fileURL, Logger log) {
    boolean canDownload = false;
    URL url;
    try {
      url = new URL(fileURL);
      HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
      int responseCode = httpConn.getResponseCode();
      if (responseCode == HttpURLConnection.HTTP_OK) {
        canDownload = true;
      }
    } catch (Exception e) {//
      log.reportException(e);
      e.printStackTrace();
    }
    return canDownload;
  }

  /**
   * Downloads a file from a URL
   *
   * @param fileURL HTTP URL of the file to be downloaded
   * @param saveFile path to where to save the file
   * @param verbose report nice things to know
   * @throws IOException
   * @return response code from HttpURLConnection
   */

  public static int downloadFile(String fileURL, String saveFile, boolean verbose,
                                 Logger log) throws IOException {
    URL url = new URL(fileURL);
    HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
    int responseCode = httpConn.getResponseCode();

    // always check HTTP response code first
    if (responseCode == HttpURLConnection.HTTP_OK) {
      new File(ext.parseDirectoryOfFile(saveFile)).mkdirs();
      String disposition = httpConn.getHeaderField("Content-Disposition");
      String contentType = httpConn.getContentType();
      int contentLength = httpConn.getContentLength();

      if (verbose) {
        log.reportTimeInfo("Content-Type = " + (contentType == null ? "unknown" : contentType));
        log.reportTimeInfo("Content-Disposition = "
                           + (disposition == null ? "unknown" : disposition));
        log.reportTimeInfo("Content-Length = " + (contentLength == -1 ? ">2GB" : contentLength));
      }
      if (contentLength < 0) {
        contentLength = Integer.MAX_VALUE;
      }

      // opens input stream from the HTTP connection
      InputStream inputStream = httpConn.getInputStream();
      String saveFilePath = saveFile;

      // opens an output stream to save into file
      FileOutputStream outputStream = new FileOutputStream(saveFilePath);

      int progressCount = 1;
      int progressVal = contentLength / 10;
      int bytesRead;
      int byteTotal = 0;
      byte[] buffer = new byte[BUFFER_SIZE];
      log.report("Downloading " + fileURL + " ", false, true);
      while ((bytesRead = inputStream.read(buffer)) != -1) {
        // consider using:
        // http://docs.oracle.com/javase/8/docs/api/javax/swing/ProgressMonitorInputStream.html
        byteTotal += bytesRead;
        if (byteTotal >= progressCount * progressVal) {
          log.report(". ", false, true);
          progressCount++;
        }
        outputStream.write(buffer, 0, bytesRead);
      }
      log.report("");

      outputStream.close();
      inputStream.close();

      log.reportTimeInfo("Downloaded " + saveFilePath);
    } else {
      log.reportError("No file to download. Server replied HTTP code: " + responseCode);
    }
    httpConn.disconnect();
    return responseCode;
  }

  public static String readFileAsHexString(String fileURL) throws IOException,
                                                           HttpDownloadException {
    // Try and open a URL connection
    URL url = new URL(fileURL);
    HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
    if (httpConn.getResponseCode() != HttpURLConnection.HTTP_OK) {
      throw new HttpDownloadException(new HTTPException(httpConn.getResponseCode()));
    }

    // Read the HTTP connection as a digest, allowing the MD5 to be computed
    try (InputStream inputStream = httpConn.getInputStream()) {
      int length = httpConn.getContentLength();
      if (length < 0) throw new HttpDownloadException("HTTP Connection (" + httpConn.toString()
                                                      + ") provided invalid content length: "
                                                      + length);
      byte[] md5 = new byte[length];
      int bytesRead = 0;

      while (bytesRead != -1) {
        int buffer = Math.min(BUFFER_SIZE, md5.length - bytesRead);
        // Read the file file
        bytesRead = inputStream.read(md5, bytesRead, buffer);
      }

      return new String(md5, ext.UTF_8).trim();
    }
  }
}

// if (disposition != null) {
// // extracts file name from header field
// int index = disposition.indexOf("filename=");
// if (index > 0) {
// fileName = disposition.substring(index + 10,
// disposition.length() - 1);
// }
// } else {
// // extracts file name from URL
// fileName = fileURL.substring(fileURL.lastIndexOf("/") + 1,
// fileURL.length());
// }
