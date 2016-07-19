package org.genvisis.link.bat;

import java.io.*;
import java.util.*;

public class osaHauser {
	public int[] maxes = {281, 253, 216, 202, 193, 180, 174, 153, 151, 168, 144, 164, 104, 127, 106, 115, 125, 126, 89, 95, 36, 47, 176};

	public osaHauser(String filename, int start, int stop, int[] newmaxes) throws IOException {
		this.maxes = newmaxes;
		new osaHauser(filename, start, stop);
	}

	public osaHauser(String filename, int start, int stop) throws IOException {
		PrintWriter writer = null;
		String chrome;

		for (int chromosome = start; chromosome<=stop; chromosome++) {
			chrome = (chromosome<10)?"0"+chromosome:""+chromosome;
			if (!(new File("chrom"+chrome)).exists()) {
				(new File("chrom"+chrome)).mkdir();
			}

			writer = new PrintWriter(new FileWriter("chrom"+chrome+"/seeds"));
			writer.println((int)(Math.random()*10000)+"\t"+(int)(Math.random()*10000)+"\t"+(int)(Math.random()*10000));
			writer.close();

			Calendar calendar = new GregorianCalendar();
			writer = new PrintWriter(new FileWriter("chrom"+chrome+"/summary"));
			writer.println((new File(".")).getAbsolutePath());
			writer.println(calendar.get(Calendar.MONTH)+"/"+calendar.get(Calendar.DAY_OF_MONTH)+"/"+calendar.get(Calendar.YEAR)+" at "+calendar.get(Calendar.HOUR_OF_DAY)+":"+calendar.get(Calendar.MINUTE));
			writer.println("npankrat");
			writer.println("morton");
			writer.println(chromosome+"   "+filename);
			writer.close();

			writer = new PrintWriter(new FileWriter("chrom"+chrome+"/osa"+chrome+".dat"));
			writer.println(filename);
			writer.println("kac.unwt.lod");
			writer.println("1   2   0.00001");
			writer.println(".");
			writer.println("0  "+maxes[chromosome-1]);
			writer.println("y y y");
			writer.println("100 10000");
			writer.println("1   1");
			writer.println("y "+filename+" Insert flippant remark here");
			writer.close();
		}

	}

	public static void main(String[] args) {
		if (args.length<1||args.length>2) {
			System.out.println("Expecting 1-2 arguments: filename and chromosome number [optional, default=all]");
			System.out.println("Note: file must contain only 2 columns - Family ID and covariate");
		} else {
			try {
				if (args.length==1) {
					new osaHauser(args[0], 1, 23);
				} else {
					new osaHauser(args[0], Integer.valueOf(args[1]).intValue(), Integer.valueOf(args[1]).intValue());
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
