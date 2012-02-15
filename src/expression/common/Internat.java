package common;

import java.net.*;
import java.io.*;

public class Internat {

	public static void printPage(String source) {
		String nextLine;
		URL url = null;
		URLConnection urlConn = null;
		InputStreamReader inStream = null;
		BufferedReader buff = null;

		try {
			// Create the URL obect that points
			// at the default file index.html
			url = new URL(source);
			urlConn = url.openConnection();
			inStream = new InputStreamReader(urlConn.getInputStream());
			buff = new BufferedReader(inStream);

			// Read and print the lines from index.html
			while (true) {
				nextLine = buff.readLine();
				if (nextLine!=null) {
					System.out.println(nextLine);
				} else {
					break;
				}
			}
		} catch (MalformedURLException e) {
			System.out.println("Please check the URL:"+e.toString());
		} catch (IOException e1) {
			System.out.println("Can't read  from the Internet: "+e1.toString());
		}
	}

	public static boolean downloadFile(String source, String output) {
		return downloadFile(source, output, true);
	}

	public static boolean downloadFile(String source, String output, boolean verbose) {
		URL remoteFile = null;
		URLConnection fileStream = null;
		DataInputStream in = null;
		// DataOutputStream out=null;
		FileOutputStream fOut = null;

		try {
			if (verbose) {
				System.out.print("Downloading "+source);
			}
			if (output.indexOf("/")>0) {
				new File(output.substring(0, output.lastIndexOf("/"))).mkdirs();
			}
			try {
				remoteFile = new URL(source);
				fileStream = remoteFile.openConnection();

				// Open the input streams for the remote file
				fOut = new FileOutputStream(output);

				// Open the output streams for saving this file on disk
				// out = new DataOutputStream(fOut);

				in = new DataInputStream(fileStream.getInputStream());
			} catch (Exception fnfe) {
				return false;
			}

			// Read the remote on save save the file
			int data;
			while ((data = in.read())!=-1) {
				fOut.write(data);
			}
			if (verbose) {
				System.out.println(" -- complete.");
			}
		} catch (Exception e) {
			// e.printStackTrace();
			return false;
		} finally {
			try {
				in.close();
				fOut.flush();
				fOut.close();
			} catch (Exception e) {
				// e.printStackTrace();
				return false;
			}

		}
		return true;
	}

	public static void main(String[] args) {

		downloadFile("http://yahoo.com", "yahoo.html");

	}
}
