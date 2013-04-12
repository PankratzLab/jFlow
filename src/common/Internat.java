package common;

import java.net.*;
import java.io.*;
import java.util.*;

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
	
//	public static void doSubmit(String url, Hashtable<String, String> data) throws Exception {
	public static String[] doSubmit(String url, Map<String, String> data, int ms_timeout) throws Exception {
		HttpURLConnection conn;
		URL siteUrl;
//		String[] keys;
		Vector<String> output;
		
		siteUrl = new URL(url);
		conn = (HttpURLConnection) siteUrl.openConnection();
		conn.setRequestMethod("POST");
		conn.setDoOutput(true);
		conn.setDoInput(true);
		
		DataOutputStream out = new DataOutputStream(conn.getOutputStream());
		
		Set<String> keys = data.keySet();
		Iterator<String> keyIter = keys.iterator();
		String content = "";
		for(int i=0; keyIter.hasNext(); i++) {
			Object key = keyIter.next();
			if(i!=0) {
				content += "&";
			}
			content += key + "=" + URLEncoder.encode(data.get(key), "UTF-8");
		}
//		System.out.println(content);
		out.writeBytes(content);
		out.flush();
		out.close();
		BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
		String line = "";
		output = new Vector<String>();
		while((line=in.readLine())!=null) {
			output.add(line);
		}
		in.close();
		
		return Array.toStringArray(output);
	}

	private static void downloadAllFilesInClipboard(String destinationDirectory) {
		String[] files;
		
		new File(destinationDirectory).mkdirs();
		files = ext.getClipboard().split("\n");
		for (int i = 0; i < files.length; i++) {
			downloadFile(files[i], destinationDirectory+ext.removeDirectoryInfo(files[i]));
		}		
	}
	
	public static void main(String[] args) {
		downloadAllFilesInClipboard("C:/ENCODE/");
//		downloadFile("http://yahoo.com", "yahoo.html");

	}
}
