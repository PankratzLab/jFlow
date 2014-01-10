package affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Vector;

import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.ext;

public class MergeChp implements Runnable{
	public static final String[][] SNP_HEADER_OPTIONS =  {{"SNP", "rsID","ProbeSetName"}};
	public static final String FILENAME_AS_ID_OPTION = "[FILENAME_ROOT]";	
	
	
	private Project proj;
	private String[] files;
	private long timeBegan;
	private int threadId;
	private String  commonSubFolderPattern;
	private String output;
	public MergeChp(Project proj, String[] files, long timeBegan) {
		this(proj, files, timeBegan , -1,"","");
	}

	public MergeChp(Project proj, String[] files,  long timeBegan, int threadId ,String commonSubFolderPattern, String output) {
		this.proj = proj;
		this.files = files;
		this.timeBegan = timeBegan;
		this.threadId = threadId;
		this.commonSubFolderPattern = commonSubFolderPattern;
		this.output = output;
	}
	
		public void run() {
			
			BufferedReader reader;
			PrintWriter writer;
			String idHeader , delimiter;
			idHeader = proj.getProperty(Project.ID_HEADER);
			delimiter = proj.getSourceFileDelimiter();
			String[] line;
			String aline;
			//common folder output by apt-genotype, present in each subdirectory of Source
			String commonSubFolder = commonSubFolderPattern;
			//check source directory
			String [] dirList = Files.listDirectories(proj.getDir(Project.SOURCE_DIRECTORY), false);
			if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
				System.err.println("Error - the Project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
				return;
			}

			int counts= 0;
			
			for (int j =0; j< files.length; j++){
				writer = Files.getAppropriateWriter(output+files[j]+".gz");
				System.out.println("merging files "+ (j+1) + " of " +  files.length);
				
				for (int i = 0; i < dirList.length; i++) {
					try{
						reader = Files.getAppropriateReader(proj.getDir(Project.SOURCE_DIRECTORY)+dirList[i]+commonSubFolder+"/"+files[j]);
						//filter comments
						do {
							line = reader.readLine().trim().split(delimiter, -1);
						} while (reader.ready()&&(ext.indexFactors(SNP_HEADER_OPTIONS, line, false, true, false, false)[0]==-1 ));
						//if its the first directory, print the header
						
						if(i ==0){
									writer.println(Array.toStr(line));
						}
						
						while (reader.ready()){
							aline = reader.readLine().trim();
							writer.println(aline);
							 counts++;
						}
						reader.close();
						
					}catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+files[j]+"\" not found in " +proj.getDir(Project.SOURCE_DIRECTORY)+dirList[i]);
						return;
					}
					catch (IOException ioe) {
						System.err.println("Error reading file \""+proj.getDir(Project.SOURCE_DIRECTORY)+dirList[i]+files[j]+"\"");
						return;
					}
				}
				writer.close();
				System.out.println("Merged file contains " + counts + " IDs");
				counts =0;
			}
		}
		
		
	public static void combineChpFiles(Project proj , int numThreads ,String commonSubFolderPattern, String commonFilename, String output )  {
		Vector<Vector<String>> fileCabinet;
		Thread[] threads;
		boolean complete;
		long timeBegan;
		//common folder output by apt-genotype, present in each subdirectory of Source
		String commonSubFolder =commonSubFolderPattern;
		//check source directory
		timeBegan = new Date().getTime();
		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
			System.err.println("Error - the Project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
			return;
		}
		

		String [] dirList = Files.listDirectories(proj.getDir(Project.SOURCE_DIRECTORY), false);
		String [] files = Files.list(proj.getDir(Project.SOURCE_DIRECTORY)+dirList[0]+commonSubFolder, commonFilename, false);
		fileCabinet = new Vector<Vector<String>>();
		for (int i = 0; i<numThreads; i++) {
			fileCabinet.add(new Vector<String>());
		}
		for (int i = 0; i<files.length; i++) {
			fileCabinet.elementAt(i%numThreads).add(files[i]);
			
		}
		System.out.println("beginning to merge " +files.length + " files");
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new MergeChp(proj, fileCabinet.elementAt(i).toArray(new String[fileCabinet.elementAt(i).size()]),  timeBegan, i ,commonSubFolder, output));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
		
		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i<numThreads; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {}
			}
		}

	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		System.out.println(numArgs);
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		String output ="";
		int numThreads = 8;
		String commonSubFolderPattern = "";
		String commonFilename = ".txt";
		String usage = "\n"+
		"affy.MergeChprequires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n";
	

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			}else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
			else if (args[i].startsWith("out=")) {
				output = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}
		//filename = "C:/workspace/Genvisis/projects/dbGaP_test_CEL.properties";
		proj = new Project(filename, false);
		try {
			combineChpFiles(proj , numThreads ,commonSubFolderPattern, commonFilename  , output );
		}
		
		catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Done Merging Files");
	}
}
