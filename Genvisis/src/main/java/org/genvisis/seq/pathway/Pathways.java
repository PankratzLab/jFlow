package seq.pathway;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import filesys.GeneData;
import filesys.GeneTrack;

public class Pathways implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final String HUMAN = "hsa";
	private static final String RS = "rs:";

	private static final String KEGG_HUMAN_PATHWAY_LIST = "http://rest.kegg.jp/list/pathway/" + HUMAN;
	private static final String KEGG_HUMAN_PATHWAY_GENE_LINK = "http://rest.kegg.jp/link/" + HUMAN + "/";
	private static final String REF_SEQ_GENE_LINK = "http://rest.genome.jp/link/refnuc/";
	private static final String PATH = "path";
	private Pathway[] pathways;
	//private Logger log;

	public Pathways(Logger log) {
		super();
		//this.log = log;
	}

	public Pathways(Pathway[] pathways, Logger log) {
		super();
		this.pathways = pathways;
		//this.log = log;
	}

	public Pathway[] getPathways() {
		return pathways;
	}

	public Pathway getPathway(String name) {
		for (int i = 0; i < pathways.length; i++) {
			if (pathways[i].getPathwayName().equals(name)) {
				return pathways[i];
			}
		}
		return null;
	}

	public Pathway[] getPathwaysFor(GeneData geneData) {
		ArrayList<Pathway> tmp = new ArrayList<Pathway>();
		for (int i = 0; i < pathways.length; i++) {
			if (pathways[i].containsGene(geneData.getGeneName())) {
				tmp.add(pathways[i]);
			}
		}
		return tmp.toArray(new Pathway[tmp.size()]);
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static Pathways load(String filename) {
		return (Pathways) Files.readSerial(filename, false, false);
	}

	private static class KeggPathwayWorker implements Callable<Pathway> {

		private String[] pathwayLine;
		private GeneTrack geneTrack;
		private Logger log;

		public KeggPathwayWorker(String[] pathwayLine, GeneTrack geneTrack, Logger log) {
			super();
			this.pathwayLine = pathwayLine;
			this.geneTrack = geneTrack;
			this.log = log;
		}

		@Override
		public Pathway call() throws Exception {
			return matchPathwayGenes(geneTrack, pathwayLine, log);
		}

	}

	private static class KeggPathwayProducer implements Producer<Pathway> {
		private GeneTrack geneTrack;
		private BufferedReader pathwayReader;
		private Logger log;

		public KeggPathwayProducer(GeneTrack geneTrack, BufferedReader pathwayReader, Logger log) {
			super();
			this.geneTrack = geneTrack;
			this.pathwayReader = pathwayReader;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			boolean ready = false;

			try {
				ready = pathwayReader.ready();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			return ready;

		}

		@Override
		public Callable<Pathway> next() {
			String[] line;
			try {
				line = pathwayReader.readLine().trim().split("\t");
				return new KeggPathwayWorker(line, geneTrack, log);

			} catch (IOException e) {
				log.reportException(e);
				e.printStackTrace();
			}
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			try {
				pathwayReader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// TODO Auto-generated method stub

		}

	}

	private void downloadKeggData(GeneTrack geneTrack, int numthreads, Logger log) {
		ArrayList<Pathway> paths = new ArrayList<Pathway>();
		BufferedReader in;
		try {
			int index = 0;
			in = new BufferedReader(new InputStreamReader(new URL(KEGG_HUMAN_PATHWAY_LIST).openStream()));
			KeggPathwayProducer producer = new KeggPathwayProducer(geneTrack, in, log);
			WorkerTrain<Pathway> train = new WorkerTrain<Pathway>(producer, numthreads, numthreads, log);
			while (train.hasNext()) {
				index++;
				if (index % 20 == 0) {
					log.reportTimeInfo("Parsed " + index + " pathways");
					// this.pathways = paths.toArray(new Pathway[paths.size()]);
					// train.shutdown();
					// return;
				}
				paths.add(train.next());
			}
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.pathways = paths.toArray(new Pathway[paths.size()]);
	}

	private static Pathway matchPathwayGenes(GeneTrack geneTrack, String[] line, Logger log) {

		HashSet<String> pathGenes = new HashSet<String>();
		String path = line[0].replaceAll(PATH + ":", "");
		String pathName = ext.replaceWithLinuxSafeCharacters(line[1], true);
		Hashtable<String, ArrayList<String>> lookup = getNCBILookupForPathway(path, log);
		GeneData[][] geneDatas = geneTrack.getGenes();
		for (String keggGene : lookup.keySet()) {
			ArrayList<String> tmp = lookup.get(keggGene);
			boolean found = false;
			for (int i = 0; i < geneDatas.length; i++) {
				if (!found) {
					for (int j = 0; j < geneDatas[i].length; j++) {
						String[] ncbis = geneDatas[i][j].getNcbiAssessionNumbers();
						for (int k = 0; k < tmp.size(); k++) {
							if (ext.indexOfStr(tmp.get(k), ncbis) >= 0) {
								found = true;
								pathGenes.add(geneDatas[i][j].getGeneName());
							}
						}
					}
				}
			}
		}
		if (pathGenes.size() != lookup.size()) {
			log.reportTimeWarning("Could not match all entries for pathway " + path + " , setting incomplete flag");
			// return new Pathway(path + ":Invalid", new GeneData[] {}, false, log);
		}
		// else {
		ArrayList<GeneData> genes = new ArrayList<GeneData>();
		for (String gene : pathGenes) {
			GeneData[] tmp = geneTrack.lookupAllGeneData(gene);
			for (int i = 0; i < tmp.length; i++) {
				genes.add(tmp[i]);
			}
		}
		return new Pathway(pathName, genes.toArray(new GeneData[genes.size()]), true, pathGenes.size() == lookup.size(), log);
	}

	private static Hashtable<String, ArrayList<String>> getNCBILookupForPathway(String pathway, Logger log) {
		ArrayList<String> keggGenes = new ArrayList<String>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new URL(KEGG_HUMAN_PATHWAY_GENE_LINK + "/" + pathway).openStream()));
			try {
				Thread.sleep(100);
			} catch (InterruptedException ie) {
			}
			while (in.ready()) {
				String[] line = in.readLine().trim().split("\t");
				String gene = line[1];
				keggGenes.add(gene);

			}
			in.close();

		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Hashtable<String, ArrayList<String>> lookup = new Hashtable<String, ArrayList<String>>();
		for (int i = 0; i < keggGenes.size(); i++) {
			lookup.put(keggGenes.get(i), new ArrayList<String>());

		}

		ArrayList<String[]> splits = Array.splitUpArray(keggGenes.toArray(new String[keggGenes.size()]), Math.round((float) keggGenes.size() / 50) + 1, log);
		for (int i = 0; i < splits.size(); i++) {

			String q = Array.toStr(splits.get(i), "+");
			String url = REF_SEQ_GENE_LINK + q;
			try {
				URLConnection uc = new URL(url).openConnection();
				uc.connect();
				BufferedReader in = new BufferedReader(new InputStreamReader(uc.getInputStream()));
				try {
					Thread.sleep(100);
				} catch (InterruptedException ie) {
				}
				while (in.ready()) {
					String[] line = in.readLine().trim().split("\t");
					String gene = line[0];
					lookup.get(gene).add(line[1].replaceAll(RS, ""));
				}
				in.close();
			} catch (MalformedURLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// log.reportTimeInfo("Using url " + url);

		}
		// if (keggGenes.size() != lookup.size()) {
		// log.reportTimeWarning("Could not look up all refseqIds for pathway " + pathway + ", " + keggGenes.size() + " genes and  found " + lookup.size());
		// for (int j = 0; j < keggGenes.size(); j++) {
		// if (!lookup.containsKey(keggGenes.get(j).replaceAll(HUMAN + ":", ""))) {
		// log.reportTimeError("Missing " + keggGenes.get(j));
		// }
		// }
		// }

		return lookup;
	}

	public static void test() {
		String geneTrack = "N:/statgen/NCBI/RefSeq_hg19.gtrack";
		GeneTrack geneTrack2 = GeneTrack.load(geneTrack, false);
		Logger log = new Logger();

		String keggSer = "D:/data/Project_Tsai_Project_021/KeggTest/kegg.ser";

		Pathways way = new Pathways(log);
		// if (!Files.exists(keggSer)) {
		way.downloadKeggData(geneTrack2, 3, log);
		way.serialize(keggSer);
		new GenomeRegions(geneTrack2, way, log).serialize(ext.addToRoot(keggSer, ".genomeRegions"));
		// } else {
		// log.reportTimeInfo("Loading pathways from " + keggSer);
		// way = load(keggSer);
		// }
	}

	public static void main(String[] args) {
		test();
	}

}