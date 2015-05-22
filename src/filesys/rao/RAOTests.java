//package filesys.rao;
//
//import java.io.FileWriter;
//import java.io.PrintWriter;
//import java.util.ArrayList;
//import java.util.Set;
//
//import cnv.var.CNVariant;
//import common.Files;
//import common.Logger;
//import filesys.rao.RAOWriter.RandomAccessProducer;
//import filesys.rao.RAOWriter.WriteComplete;
//
//public class RAOTests {
//
//	public static void test(String testFile) {
//		Logger log = new Logger(testFile + ".log");
//		CNVariant[] segs = new CNVariant[] { new CNVariant("HFS", "DS", (byte) 1, 32, 32, 3, 33, 10, -1) };
//		ArrayList<CNVariant> daSegs = new ArrayList<CNVariant>();
//		for (int i = 0; i < 100000; i++) {
//			daSegs.add(new CNVariant("HFS", "DS", (byte) 1, i, i + 100, 3, 33, 10, -1));
//		}
//		segs = daSegs.toArray(new CNVariant[daSegs.size()]);
//		final CNVariant[] segsdump = segs;
//
//		RandomAccessProducer producer = new RandomAccessProducer() {
//			int index = 0;
//
//			@Override
//			public void remove() {
//				// TODO Auto-generated method stub
//
//			}
//
//			@Override
//			public RAObject next() {
//				final CNVariant Aseg = segsdump[index];
//				index++;
//				return Aseg;
//			}
//
//			@Override
//			public boolean hasNext() {
//				return index < segsdump.length;
//
//			}
//		};
//
//		RAOWriter writer = new RAOWriter(testFile, producer, log);
//		log.reportTimeInfo("Beginning ser write");
//		long time = System.currentTimeMillis();
//		WriteComplete writeComplete = writer.writeToFile();
//		RAOIndex raoIndex = RAOIndex.load(writeComplete.getIndexFile(), log);
//		Set<String> s = raoIndex.getIndex().keySet();
//		RAOReader<CNVariant> reader = new RAOReader<CNVariant>(writeComplete.getDataFile(), writeComplete.getIndexFile(), log);
//		for (String as : s) {
//			System.out.println(as + "\t" + reader.loadPosition(raoIndex.getIndex().get(as).get(1000)).getUCSClocation());
//		}
//
//		log.reportTimeInfo("Ending ser write");
//		log.reportTimeElapsed(time);
//
//		try {
//			PrintWriter awriter = new PrintWriter(new FileWriter(testFile + ".txt"));
//			log.reportTimeInfo("Beginning normal write");
//
//			time = System.currentTimeMillis();
//			for (int i = 0; i < segsdump.length; i++) {
//				awriter.println(segsdump[i].toPlinkFormat());
//			}
//			log.reportTimeInfo("Ending normal write");
//			log.reportTimeElapsed(time);
//			awriter.close();
//		} catch (Exception e) {
//			log.reportError("Error writing to " + testFile + ".txt");
//			log.reportException(e);
//		}
//
//		log.reportTimeInfo("Beginning normal ser");
//
//		time = System.currentTimeMillis();
//		Files.writeSerial(segsdump, testFile + ".regSer", true);
//
//		log.reportTimeInfo("Ending normal ser");
//		log.reportTimeElapsed(time);
//	}
//
//	public static void main(String[] args) {
//		String testFile = "D:/data/Project_Tsai_21_25_26_spector/testRanAccs/ran2.rao";
//		test(testFile);
//	}
//
//}
