package common;

import java.io.*;
import java.util.*;
import java.util.zip.*;

public class Zip {

	public static void zip(String filename) {
		zip(new String[] {filename}, filename+".zip");
	}

	public static void zip(String[] filenames, String zipfile) {
		byte[] buf = new byte[1024];

		try {
			ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipfile));

			for (int i = 0; i<filenames.length; i++) {
				FileInputStream in = new FileInputStream(filenames[i]);

				out.putNextEntry(new ZipEntry(filenames[i]));
				int len;
				while ((len = in.read(buf))>0) {
					out.write(buf, 0, len);
				}

				out.closeEntry();
				in.close();
			}

			out.close();
		} catch (IOException e) {
			System.err.println("Error creating zipfile '"+zipfile+"'");
			e.printStackTrace();
		}

	}

	public static void zipDirectory(String dir, String outfile) {
		try {
			ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(outfile));
			zipDir(dir, zos);
			zos.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void zipDir(String dir2zip, ZipOutputStream zos) {
		try {
			File zipDir = new File(dir2zip);
			String[] dirList = zipDir.list();
			byte[] readBuffer = new byte[2156];
			int bytesIn = 0;

			for (int i = 0; i<dirList.length; i++) {
				File f = new File(zipDir, dirList[i]);
				if (f.isDirectory()) {
					String filePath = f.getPath();
					zipDir(filePath, zos);
				} else {
					FileInputStream fis = new FileInputStream(f);
					String name = f.getPath();
					if (name.startsWith(".")) {
						name = name.substring(2);
					}
					ZipEntry anEntry = new ZipEntry(name);
					zos.putNextEntry(anEntry);
					while ((bytesIn = fis.read(readBuffer))!=-1) {
						zos.write(readBuffer, 0, bytesIn);
					}
					fis.close();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void unzipFile(String zipfile, String destination) {
		if (zipfile.endsWith(".gz")) {
			gunzip(zipfile, destination);
		} else {
			unzip(zipfile, destination);
		}
	}
	
	private static void unzip(String zipfile, String destination) {
		try {
			ZipInputStream in = new ZipInputStream(new FileInputStream(zipfile));

			if (!destination.endsWith("/")&&!destination.endsWith("\\")) {
				destination += "/";
			}

			String[] files = list(zipfile);
			for (int i = 0; i<files.length; i++) {
				ZipEntry entry = in.getNextEntry();
				String filename = destination+entry.getName();

				if (!filename.endsWith("/")&&!filename.endsWith("\\")) {
					int index = Math.max(filename.lastIndexOf("/"), filename.lastIndexOf("\\"));

					if (index>=0) {
						new File(filename.substring(0, index)).mkdirs();
					}

					OutputStream out = new FileOutputStream(destination+entry.getName());

					byte[] buf = new byte[1024];
					int len;
					while ((len = in.read(buf))>0) {
						out.write(buf, 0, len);
					}

					out.close();
				}

			}

			in.close();
		} catch (IOException e) {
			System.err.println("Error unzipping '"+zipfile+"'");
			e.printStackTrace();
		}

	}

	public static void gzip(String filename) {
		gzip(filename, filename+".gz");
	}
	
	public static void gzip(String filename, String outputFilename) {
		GZIPOutputStream out;
		FileInputStream in;
		byte[] buf = new byte[1024];
		int len;

		try {
			out = new GZIPOutputStream(new FileOutputStream(outputFilename));
			in = new FileInputStream(filename);

			while ((len = in.read(buf))>0) {
				out.write(buf, 0, len);
			}
			in.close();
			out.close();
		} catch (IOException e) {
			System.err.println("Error creating zipfile '"+outputFilename+"'");
			e.printStackTrace();
		}
	}

	private static void gunzip(String zipfile, String destinationDirectory) {
		try {
			GZIPInputStream in = new GZIPInputStream(new FileInputStream(zipfile));

			if (!destinationDirectory.endsWith("/")&&!destinationDirectory.endsWith("\\")) {
				destinationDirectory += "/";
			}

			OutputStream out = new FileOutputStream(destinationDirectory+ext.rootOf(zipfile, true));

			byte[] buf = new byte[1024];
			int len;
			while ((len = in.read(buf))>0) {
				out.write(buf, 0, len);
			}

			out.close();

			in.close();
		} catch (IOException e) {
			System.err.println("Error unzipping '"+zipfile+"'");
			e.printStackTrace();
		}

	}
	
	public static void gzipAllTextFilesInDirectory(String dirin, String dirout, Logger log) {
		String[] files;
		
		
		if (dirin == null) {
			log.reportError("Error - dir cannot be null");
			return;
		}

		if (!Files.exists(dirin)) {
			log.reportError("Error - invalid directory: "+dirin);
			return;
		}
		
		dirin = ext.verifyDirFormat(dirin);
		if (dirout == null) {
			dirout = dirin+"compressed/";
		}
		
		new File(dirout).mkdirs();
		files = Files.list(dirin, null, false);
		for (int i = 0; i < files.length; i++) {
			log.report("compressing file #"+(i+1)+" of "+files.length+" ("+files[i]+")");
			gzip(dirin+files[i], dirout+files[i]+".gz");
		}
	}
	
	public static void unzipParticularFile(String target, String zipfile, String destination) {
		unzipParticularFiles(new String[] {target}, zipfile, destination);
	}

	public static void unzipParticularFiles(String[] targets, String zipfile, String destination) {
		try {
			ZipInputStream in = new ZipInputStream(new FileInputStream(zipfile));

			if (!destination.endsWith("/")&&!destination.endsWith("\\")) {
				destination += "/";
			}

			String[] files = list(zipfile);
			int numErrors = 0;
			for (int i = 0; i<targets.length; i++) {
				// targets[i] = ext.replaceAllWith(targets[i], "\\", "/");
				if (ext.indexOfStr(targets[i], files)==-1) {
					System.err.println("Error - '"+targets[i]+"' was not found in '"+zipfile+"'");
				}
			}
			if (numErrors>0) {
				System.exit(1);
			}
			for (int i = 0; i<files.length; i++) {
				ZipEntry entry = in.getNextEntry();
				String filename = destination+entry.getName();

				if (ext.indexOfStr(entry.getName(), targets)!=-1&&!filename.endsWith("/")&&!filename.endsWith("\\")) {
					int index = Math.max(filename.lastIndexOf("/"), filename.lastIndexOf("\\"));

					if (index>=0) {
						new File(filename.substring(0, index)).mkdirs();
					}

					OutputStream out = new FileOutputStream(destination+entry.getName());

					byte[] buf = new byte[1024];
					int len;
					while ((len = in.read(buf))>0) {
						out.write(buf, 0, len);
					}

					out.close();
				}

			}

			in.close();
		} catch (IOException e) {
			System.err.println("Error unzipping '"+zipfile+"'");
			e.printStackTrace();
		}

	}

	@SuppressWarnings("unchecked")
	public static String[] list(String filename) {
		Vector<String> v = new Vector<String>();

		try {
			ZipFile zf = new ZipFile(filename);
			for (Enumeration<ZipEntry> entries = (Enumeration<ZipEntry>)zf.entries(); entries.hasMoreElements();) {
				v.add((entries.nextElement()).getName());
			}
			zf.close();
		} catch (IOException e) {
			System.err.println("Error listing contents of zipfile '"+filename+"'");
			e.printStackTrace();
		}

		return Array.toStringArray(v);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String dirin = "./";
		String dirout = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + 
		"common.Zip requires 0-1 arguments\n" + 
		"   (1) name of file to be gzipped (i.e. file=in.csv (not the default))\n" + 
		" OR:\n" + 
		"   (1) directory containing files to be gzipped (i.e. dirin=" + dirin + " (default))\n" + 
		"   (2) directory where the compressed files should be written to (i.e. dirout=[dirin]/compressed/ (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dirin=")) {
				dirin = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dirout=")) {
				dirin = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			
//			dirin = "I:/GEDI/02_After_reclustering/00src/";
//			dirout = "C:/GEDI_compressed/";
			
//			dirin = "D:/data/GEDI/penn_data/";
			
//			dirin = "D:/data/genomestudio/GEDI_Exome/FinalReportsChargeTop/";
//			dirout = "C:/GEDI_exomeChargeTop_compressed/";

			dirin = "C:/LoganReports/";
			dirout = "I:/Logan/";
			
			log = new Logger(logfile);
			if (filename != null) {
				gzip(filename);
			} else if (dirin != null) {
				gzipAllTextFilesInDirectory(dirin, dirout, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
